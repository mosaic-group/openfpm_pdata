#ifndef CARTESIAN_GRAPH_UNIT_TEST_HPP
#define CARTESIAN_GRAPH_UNIT_TEST_HPP

#include "Graph/CartesianGraphFactory.hpp"
#include "Graph/map_graph.hpp"

#define GS_SIZE 8

/*!
 *
 * Test node
 *
 */

struct node_cp
{
	//! The node contain 3 unsigned long integer for comunication computation and memory
	typedef boost::fusion::vector<size_t,size_t,size_t> type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! communication property id in boost::fusion::vector
	static const unsigned int communication = 0;
	//! computation property id in boost::fusion::vector
	static const unsigned int computation = 1;
	//! memory property id in boost::fusion::vector
	static const unsigned int memory = 2;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 3;
};

const std::string node_cp::attributes::name[] = {"communication","computation","memory"};

BOOST_AUTO_TEST_SUITE( CartesianGraphFactory_test )

BOOST_AUTO_TEST_CASE( CartesianGraphFactory_use_np)
{
	typedef node_cp node;

	CartesianGraphFactory<3,Graph_CSR<Point_test<float>,Point_test<float>>> g_factory;

	// Cartesian grid
	size_t sz[3] = {GS_SIZE,GS_SIZE,GS_SIZE};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions, non periodic
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	Graph_CSR<Point_test<float>,Point_test<float>> g = g_factory.construct<node::communication,NO_VERTEX_ID,float,2>(sz,box,bc);

	// check that the number of vertex are equal to GS_SIZE^3
	BOOST_REQUIRE_EQUAL(g.getNVertex(),(size_t)GS_SIZE*GS_SIZE*GS_SIZE);

	// check that the number of vertex are equal to GS_SIZE^3
	BOOST_REQUIRE_EQUAL(g.getNEdge(),(size_t)3*8+4*(GS_SIZE-2)*12+6*(GS_SIZE-2)*(GS_SIZE-2)*5+(GS_SIZE-2)*(GS_SIZE-2)*(GS_SIZE-2)*6);
}

BOOST_AUTO_TEST_CASE( CartesianGraphFactory_use_p)
{
	typedef node_cp node;

	CartesianGraphFactory<3,Graph_CSR<Point_test<float>,Point_test<float>>> g_factory;

	// Cartesian grid
	size_t sz[3] = {GS_SIZE,GS_SIZE,GS_SIZE};

	// Box
	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions, non periodic
	size_t bc[] = {PERIODIC,PERIODIC,PERIODIC};

	Graph_CSR<Point_test<float>,Point_test<float>> g = g_factory.construct<node::communication,NO_VERTEX_ID,float,2>(sz,box,bc);

	// check that the number of vertex are equal to GS_SIZE^3
	BOOST_REQUIRE_EQUAL(g.getNVertex(),(size_t)GS_SIZE*GS_SIZE*GS_SIZE);

	// check that the number of vertex are equal to GS_SIZE^3
	BOOST_REQUIRE_EQUAL(g.getNEdge(),(size_t)(GS_SIZE)*(GS_SIZE)*(GS_SIZE)*6);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
