#ifndef DIST_MAP_GRAPH_UNIT_TEST_HPP
#define DIST_MAP_GRAPH_UNIT_TEST_HPP

#include "Graph/DistGraphFactory.hpp"
#include "Graph/dist_map_graph.hpp"
#include "Packer.hpp"
#include "Unpacker.hpp"

#define GS_SIZE 4

struct vx
{
	typedef boost::fusion::vector<float[3], size_t, size_t, size_t, size_t, size_t> type;

	typedef typename memory_traits_inte<type>::type memory_int;
	typedef typename memory_traits_lin<type>::type memory_lin;

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
	//! id property id in boost::fusion::vector
	static const unsigned int id = 1;
	//! global id property id in boost::fusion::vector
	static const unsigned int global_id = 2;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 3;

	vx()
	{

	}

	inline vx(const vx & p)
	{
		boost::fusion::at_c<0>(data)[0] = boost::fusion::at_c<0>(p.data)[0];
		boost::fusion::at_c<0>(data)[1] = boost::fusion::at_c<0>(p.data)[1];
		boost::fusion::at_c<0>(data)[2] = boost::fusion::at_c<0>(p.data)[2];
		boost::fusion::at_c<1>(data) = boost::fusion::at_c<1>(p.data);
		boost::fusion::at_c<2>(data) = boost::fusion::at_c<2>(p.data);
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
		boost::fusion::at_c<1>(data) = p.template get<1>();
		boost::fusion::at_c<2>(data) = p.template get<2>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

const std::string vx::attributes::name[] = { "x", "id", "global_id" };

struct ed
{
	typedef boost::fusion::vector<size_t, size_t> type;

	//! Attributes name
	struct attributes
	{
		static const std::string name[];
	};

	//! The data
	type data;

	//! srcgid property id in boost::fusion::vector
	static const unsigned int srcgid = 0;
	//! dstgid property id in boost::fusion::vector
	static const unsigned int dstgid = 1;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 2;

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
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};

const std::string ed::attributes::name[] = { "srcgid", "dstgid" };

BOOST_AUTO_TEST_SUITE (dist_map_graph_test)

BOOST_AUTO_TEST_CASE( dist_map_graph_use_4p)
{

	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

	//! Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! Cartesian grid
	size_t sz[2] = { GS_SIZE, GS_SIZE };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph factory
	DistGraphFactory<2, DistGraph_CSR<vx, ed>> g_factory;

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd = g_factory.template construct<NO_EDGE, vx::id, vx::global_id, ed::srcgid, ed::dstgid, float, 2 - 1, 0, 1, 2>(sz, box);

	//! [Request some vertices given global ids]

	gd.deleteGhosts();

	gd.reqVertex(13);
	gd.reqVertex(1);
	gd.reqVertex(14);
	gd.sync();

	BOOST_REQUIRE_EQUAL(gd.getVertex(13).template get<vx::global_id>(), 13);
	BOOST_REQUIRE_EQUAL(gd.getVertex(1).template get<vx::global_id>(), 1);
	BOOST_REQUIRE_EQUAL(gd.getVertex(14).template get<vx::global_id>(), 14);

	gd.reqVertex(15);
	gd.reqVertex(2);
	gd.reqVertex(10);
	gd.sync();

	BOOST_REQUIRE_EQUAL(gd.getVertex(15).template get<vx::global_id>(), 15);
	BOOST_REQUIRE_EQUAL(gd.getVertex(2).template get<vx::global_id>(), 2);
	BOOST_REQUIRE_EQUAL(gd.getVertex(10).template get<vx::global_id>(), 10);

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

	/*
	 for (size_t i = 0; i < gd.getNVertex(); i++)
	 {
	 for (size_t s = 0; s < gd.getNChilds(i); ++s)
	 {
	 size_t vid = gd.getChildEdge(i, s).template get<ed::dstgid>();
	 std::cout << vid << " ";
	 }
	 }
	 std::cout << "\n";*/

	if(vcl.getProcessUnitID() == 0)
	{
		BOOST_REQUIRE_EQUAL(gd.vertex(0).template get<vx::id>(), 0);
		BOOST_REQUIRE_EQUAL(gd.vertex(1).template get<vx::id>(), 1);
		BOOST_REQUIRE_EQUAL(gd.vertex(2).template get<vx::id>(), 2);
		BOOST_REQUIRE_EQUAL(gd.vertex(3).template get<vx::id>(), 3);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(0), 3);
		BOOST_REQUIRE_EQUAL(gd.getChild(0,0), 1);
		BOOST_REQUIRE_EQUAL(gd.getChild(0,1), 14);
		BOOST_REQUIRE_EQUAL(gd.getChild(0,2), 4);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(1), 4);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,0), 2);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,1), 0);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,2), 15);
		BOOST_REQUIRE_EQUAL(gd.getChild(1,3), 5);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(2), 4);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,0), 3);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,1), 1);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,2), 8);
		BOOST_REQUIRE_EQUAL(gd.getChild(2,3), 6);
		BOOST_REQUIRE_EQUAL(gd.getNChilds(3), 3);
		BOOST_REQUIRE_EQUAL(gd.getChild(3,0), 2);
		BOOST_REQUIRE_EQUAL(gd.getChild(3,1), 9);
		BOOST_REQUIRE_EQUAL(gd.getChild(3,2), 7);
	}

	if(vcl.getProcessUnitID() == 1)
	{
		BOOST_REQUIRE_EQUAL(gd.vertex(0).template get<vx::id>(), 4);
		BOOST_REQUIRE_EQUAL(gd.vertex(1).template get<vx::id>(), 5);
		BOOST_REQUIRE_EQUAL(gd.vertex(2).template get<vx::id>(), 6);
		BOOST_REQUIRE_EQUAL(gd.vertex(3).template get<vx::id>(), 7);
	}

	if(vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(gd.vertex(0).template get<vx::id>(), 8);
		BOOST_REQUIRE_EQUAL(gd.vertex(1).template get<vx::id>(), 9);
		BOOST_REQUIRE_EQUAL(gd.vertex(2).template get<vx::id>(), 10);
		BOOST_REQUIRE_EQUAL(gd.vertex(3).template get<vx::id>(), 11);
		BOOST_REQUIRE_EQUAL(gd.vertex(4).template get<vx::id>(), 12);
		BOOST_REQUIRE_EQUAL(gd.vertex(5).template get<vx::id>(), 13);
	}

	if(vcl.getProcessUnitID() == 3)
	{
		BOOST_REQUIRE_EQUAL(gd.vertex(0).template get<vx::id>(), 14);
		BOOST_REQUIRE_EQUAL(gd.vertex(1).template get<vx::id>(), 15);
	}
}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_4p_redistribution)
{
	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

	//! Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! Cartesian grid
	size_t sz[2] = { 4, 4 };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph factory
	DistGraphFactory<2, DistGraph_CSR<vx, ed>> g_factory;

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd = g_factory.template construct<NO_EDGE, vx::id, vx::global_id, ed::srcgid, ed::dstgid, float, 2 - 1, 0, 1, 2>(sz, box);

	for(size_t i=0; i< gd.getNVertex(); i++)
	gd.vertex(i).template get<vx::x>()[2] = 0;

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

	if (vcl.getProcessUnitID() == 0)
	{
		gd.deleteGhosts();

		for (size_t i = 0; i < gd.getTotNVertex(); ++i)
		{
			gd.reqVertex(i);
		}
	}

	gd.sync();

	if (vcl.getProcessUnitID() == 0)
	{
		VTKWriter<DistGraph_CSR<vx, ed>, GRAPH> gv2(gd);
		gv2.write("dist_graph_0.vtk");

		bool test = compare("dist_graph_0.vtk","dist_graph_0_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);

		gd.deleteGhosts();
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

	if (vcl.getProcessUnitID() == 0)
	{
		gd.deleteGhosts();

		for (size_t i = 0; i < gd.getTotNVertex(); ++i)
		{
			gd.reqVertex(i);
		}
	}

	gd.sync();

	if (vcl.getProcessUnitID() == 0)
	{
		VTKWriter<DistGraph_CSR<vx, ed>, DIST_GRAPH> gv2(gd);
		gv2.write("dist_graph_1.vtk");

		bool test = compare("dist_graph_1.vtk","dist_graph_1_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);

		gd.deleteGhosts();
	}

}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_4p_free_add)
{
	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

	//! Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! Cartesian grid
	size_t sz[2] = { 4, 4 };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd;

	for (size_t i = 0; i < 4; ++i)
	{
		vx v;
		gd.addVertex(v);
	}

	gd.initProperties();

	if(vcl.getProcessUnitID() == 0)
	{
		for (size_t i = 0; i < 4; ++i)
		{
			BOOST_REQUIRE_EQUAL(gd.vertex(i).template get<vx::id>(), i);
		}
	}

	if(vcl.getProcessUnitID() == 2)
	{
		for (size_t i = 8, j = 0; i < 12 && j < gd.getNVertex(); ++i, ++j)
		{
			BOOST_REQUIRE_EQUAL(gd.vertex(j).template get<vx::id>(), i);
		}
	}

	for(size_t i=0; i< gd.getNVertex(); i++)
		gd.vertex(i).template get<vx::x>()[2] = 0;

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

	if (vcl.getProcessUnitID() == 0)
	{
		gd.deleteGhosts();

		for (size_t i = 0; i < gd.getTotNVertex(); ++i)
		{
			gd.reqVertex(i);
		}
	}

	gd.sync();

	if (vcl.getProcessUnitID() == 0)
	{
		VTKWriter<DistGraph_CSR<vx, ed>, GRAPH> gv2(gd);
		gv2.write("dist_graph_0_free.vtk");

		bool test = compare("dist_graph_0_free.vtk","dist_graph_0_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);

		gd.deleteGhosts();
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

	if (vcl.getProcessUnitID() == 0)
	{
		gd.deleteGhosts();

		for (size_t i = 0; i < gd.getTotNVertex(); ++i)
		{
			gd.reqVertex(i);
		}
	}

	gd.sync();

	if (vcl.getProcessUnitID() == 0)
	{
		VTKWriter<DistGraph_CSR<vx, ed>, DIST_GRAPH> gv2(gd);
		gv2.write("dist_graph_1_free.vtk");

		bool test = compare("dist_graph_1_free.vtk","dist_graph_1_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);

		gd.deleteGhosts();
	}

}

BOOST_AUTO_TEST_SUITE_END()

#endif
