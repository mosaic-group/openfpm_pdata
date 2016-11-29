#ifndef DIST_MAP_GRAPH_UNIT_TEST_HPP
#define DIST_MAP_GRAPH_UNIT_TEST_HPP

#include "Graph/DistGraphFactory.hpp"
#include "Graph/dist_map_graph.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"

#define DGRAPH_GS_SIZE 4

struct vx
{
	typedef boost::fusion::vector<float[3]> type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! type of the positional field
	typedef float s_type;

	//! The data
	type data;

	//! x property id in boost::fusion::vector
	static const unsigned int x = 0;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 1;

	vx()
	{

	}

	inline vx(const vx & p)
	{
		boost::fusion::at_c<0>(data)[0] = boost::fusion::at_c<0>(p.data)[0];
		boost::fusion::at_c<0>(data)[1] = boost::fusion::at_c<0>(p.data)[1];
		boost::fusion::at_c<0>(data)[2] = boost::fusion::at_c<0>(p.data)[2];
	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int id> inline auto get() const -> const decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int dim, typename Mem> inline vx(const encapc<dim, vx, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline vx & operator=(const encapc<dim, vx, Mem> & p)
	{
		boost::fusion::at_c<0>(data)[0] = p.template get<0>()[0];
		boost::fusion::at_c<0>(data)[1] = p.template get<0>()[1];
		boost::fusion::at_c<0>(data)[2] = p.template get<0>()[2];

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

const std::string vx::attributes::name[] = { "x" };

struct ed
{
	typedef boost::fusion::vector<size_t> type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! srcgid property id in boost::fusion::vector
	static const unsigned int prop = 0;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 1;

	ed()
	{

	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int dim, typename Mem> inline ed(const encapc<dim, ed, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline ed & operator=(const encapc<dim, ed, Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.get<0>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

const std::string ed::attributes::name[] = { "prop" };

BOOST_AUTO_TEST_SUITE (dist_map_graph_test)

BOOST_AUTO_TEST_CASE( dist_map_graph_use)
{

	//! Vcluster
	Vcluster & vcl = create_vcluster();

	if(vcl.getProcessingUnits() != 4)
		return;

	//! Cartesian grid
	size_t sz[2] = { DGRAPH_GS_SIZE, DGRAPH_GS_SIZE };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph factory
	DistGraphFactory<2, DistGraph_CSR<vx, ed>> g_factory;

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd = g_factory.construct<NO_EDGE, float, 2 - 1, 0, 1, 2>(sz, box);

	//! [Request some vertices given global ids]

	gd.reqVertex(13);
	gd.reqVertex(1);
	gd.reqVertex(14);
	gd.sync();

	gd.reqVertex(15);
	gd.reqVertex(2);
	gd.reqVertex(10);
	gd.sync();

	gd.deleteGhosts();

	//! [Exchange n vertices and edges packed]

	if(vcl.getProcessUnitID() == 0)
	{
		for(size_t i = 0; i < 4; i++)
		gd.q_move(i, 1);
	}

	if(vcl.getProcessUnitID() == 1)
	{
		for(size_t i = 0; i < 4; i++)
		gd.q_move(i, 0);
	}

	if(vcl.getProcessUnitID() == 2)
	{
		for(size_t i = 0; i < 2; i++)
		gd.q_move(i, 3);
	}

	if(vcl.getProcessUnitID() == 3)
	{
		for(size_t i = 0; i < 4; i++)
		gd.q_move(i, 2);
	}

	//! Redistribute
	gd.redistribute();

	if(vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 0ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 1ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(2), 2ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(3), 3ul);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(0), 3ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(0,0), 1ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(0,1), 14ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(0,2), 4ul);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(1), 4ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,0), 2ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,1), 0ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,2), 15ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,3), 5ul);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(2), 4ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,0), 3ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,1), 1ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,2), 8ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,3), 6ul);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(3), 3ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(3,0), 2ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(3,1), 9ul);
		BOOST_REQUIRE_EQUAL(gd.getChild(3,2), 7ul);
	}

	if(vcl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 4ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 5ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(2), 6ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(3), 7ul);
	}

	if(vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 8ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 9ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(2), 10ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(3), 11ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(4), 12ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(5), 13ul);
	}

	if(vcl.getProcessUnitID() == 3)
	{
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 14ul);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 15ul);
	}
}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_redistribution)
{
	//! Vcluster
	Vcluster & vcl = create_vcluster();

	if(vcl.getProcessingUnits() != 4)
		return;

	//! Cartesian grid
	size_t sz[2] = { 4, 4 };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph factory
	DistGraphFactory<2, DistGraph_CSR<vx, ed>> g_factory;

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd = g_factory.construct<NO_EDGE, float, 2 - 1, 0, 1, 2>(sz, box);

	for(size_t i=0; i< gd.getNVertex(); i++)
	gd.vertex(i).get<vx::x>()[2] = 0;

	if (vcl.getProcessUnitID() == 0)
	{
		gd.q_move(0,1);
		gd.q_move(1,1);
	}

	if (vcl.getProcessUnitID() == 1)
	{
		gd.q_move(2,0);
		gd.q_move(3,0);
	}

	if (vcl.getProcessUnitID() == 2)
	{
		gd.q_move(0,3);
		gd.q_move(1,3);
	}

	if (vcl.getProcessUnitID() == 3)
	{
		gd.q_move(2,2);
		gd.q_move(3,2);
	}

	gd.redistribute();

	VTKWriter<DistGraph_CSR<vx, ed>, DIST_GRAPH> gv2(gd);
	gv2.write("dist_graph_redistribution_0.vtk");

	if(vcl.getProcessUnitID() == 0)
	{
		bool test = compare("dist_graph_redistribution_0.vtk", "src/Graph/test_data/dist_graph_redistribution_0_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}

	if (vcl.getProcessUnitID() == 2)
	{
		gd.q_move(0,0);
		gd.q_move(1,0);
	}

	if (vcl.getProcessUnitID() == 3)
	{
		gd.q_move(0,0);
		gd.q_move(2,0);
	}

	gd.redistribute();

	gv2.write("dist_graph_redistribution_1.vtk");

	if(vcl.getProcessUnitID() == 0)
	{
		bool test = compare("dist_graph_redistribution_1.vtk","src/Graph/test_data/dist_graph_redistribution_1_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}

}


BOOST_AUTO_TEST_CASE( dist_map_graph_use_free_add)
{
	// Vcluster
	Vcluster & vcl = create_vcluster();

	if(vcl.getProcessingUnits() != 4)
		return;

	// [create graph adding freely the vertices and the edges ]

	// Distributed graph
	DistGraph_CSR<vx, ed> gd;

	// Add vertices
	for (size_t i = 0; i < 4; ++i)
	{
		vx v;
		v.get<vx::x>()[0] = vcl.getProcessUnitID();
		v.get<vx::x>()[1] = i;
		v.get<vx::x>()[2] = 0;
		size_t gid = vcl.getProcessUnitID()*4 + i;
		gd.add_vertex(v, gid);
	}

	// This method must be called after adding vertices
	gd.init();

	// Add edges
	if(vcl.getProcessUnitID()==0)
	{
		gd.add_edge(0,1);
		gd.add_edge(1,2);
		gd.add_edge(2,3);
		gd.add_edge(0,4);
		gd.add_edge(1,5);
		gd.add_edge(2,6);
		gd.add_edge(3,7);
	}

	if(vcl.getProcessUnitID()==1)
	{
		gd.add_edge(4,5);
		gd.add_edge(5,6);
		gd.add_edge(6,7);
		gd.add_edge(4,8);
		gd.add_edge(5,9);
		gd.add_edge(6,10);
		gd.add_edge(7,11);
	}

	if(vcl.getProcessUnitID()==2)
	{
		gd.add_edge(8,9);
		gd.add_edge(9,10);
		gd.add_edge(10,11);
		gd.add_edge(8,12);
		gd.add_edge(9,13);
		gd.add_edge(10,14);
		gd.add_edge(11,15);
	}

	if(vcl.getProcessUnitID()==3)
	{
		gd.add_edge(12,13);
		gd.add_edge(13,14);
		gd.add_edge(14,15);
	}

	gd.syncEdge();

	//! [create graph adding freely the vertices and the edges ]

	if(vcl.getProcessUnitID() == 0)
	{
		for (size_t i = 0; i < 4; ++i)
		{
			BOOST_REQUIRE_EQUAL(gd.getVertexId(i), i);
		}
	}

	if(vcl.getProcessUnitID() == 2)
	{
		for (size_t i = 8, j = 0; i < 12 && j < gd.getNVertex(); ++i, ++j)
		{
			BOOST_REQUIRE_EQUAL(gd.getVertexId(j), i);
		}
	}

	if(vcl.getProcessUnitID() == 0)
		gd.reqVertex(5);

	gd.sync();

	if(vcl.getProcessUnitID() == 0)
		BOOST_REQUIRE_EQUAL(gd.getVertexId(4), 5ul);

	gd.deleteGhosts();

	if (vcl.getProcessUnitID() == 0)
	{
		gd.q_move(0,1);
		gd.q_move(1,1);
	}

	if (vcl.getProcessUnitID() == 1)
	{
		gd.q_move(2,0);
		gd.q_move(3,0);
	}

	if (vcl.getProcessUnitID() == 2)
	{
		gd.q_move(0,3);
		gd.q_move(1,3);
	}

	if (vcl.getProcessUnitID() == 3)
	{
		gd.q_move(2,2);
		gd.q_move(3,2);
	}

	gd.redistribute();

	VTKWriter<DistGraph_CSR<vx, ed>, DIST_GRAPH> gv2(gd);
	gv2.write("dist_graph_free_0.vtk");

	if(vcl.getProcessUnitID() == 0)
	{
		bool test = compare("dist_graph_free_0.vtk", "src/Graph/test_data/dist_graph_free_0_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}

	if (vcl.getProcessUnitID() == 2)
	{
		gd.q_move(0,0);
		gd.q_move(1,0);
	}

	if (vcl.getProcessUnitID() == 3)
	{
		gd.q_move(0,0);
		gd.q_move(2,0);
	}

	gd.redistribute();

	gv2.write("dist_graph_free_1.vtk");

	if(vcl.getProcessUnitID() == 0)
	{
		bool test = compare("dist_graph_free_1.vtk", "src/Graph/test_data/dist_graph_free_1_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}
}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_multi_free_add)
{
	// Vcluster
	Vcluster & vcl = create_vcluster();

	if(vcl.getProcessingUnits() != 4)
		return;

	// Distributed graph
	DistGraph_CSR<vx, ed> gd;

	// Add vertices
	if(vcl.getProcessUnitID()==0)
	{
		vx v;
		gd.add_vertex(v, 0);
		gd.add_vertex(v, 1);
		gd.add_vertex(v, 2);
		gd.add_vertex(v, 3);
	}

	if(vcl.getProcessUnitID()==1)
	{
		vx v;
		gd.add_vertex(v, 4);
		gd.add_vertex(v, 5);
		gd.add_vertex(v, 6);
		gd.add_vertex(v, 7);
	}

	// This method must be called ALWAYS after adding vertices
	gd.init();

	if(vcl.getProcessUnitID()==2)
	{
		vx v;
		gd.add_vertex(v, 8);
		gd.add_vertex(v, 9);
		gd.add_vertex(v, 10);
		gd.add_vertex(v, 11);
	}

	if(vcl.getProcessUnitID()==3)
	{
		vx v;
		gd.add_vertex(v, 12);
		gd.add_vertex(v, 13);
		gd.add_vertex(v, 14);
		gd.add_vertex(v, 15);
	}

	// This method must be called ALWAYS after adding vertices
	gd.init();

	gd.reqVertex(15);

	gd.sync();

	BOOST_REQUIRE_EQUAL(gd.getVertexId(gd.getNVertex()-1), 15ul);
}

BOOST_AUTO_TEST_SUITE_END()

#endif
