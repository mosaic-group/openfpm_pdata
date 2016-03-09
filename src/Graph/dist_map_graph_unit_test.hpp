#ifndef DIST_MAP_GRAPH_UNIT_TEST_HPP
#define DIST_MAP_GRAPH_UNIT_TEST_HPP

#include "Graph/DistGraphFactory.hpp"
#include "Graph/dist_map_graph.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"

#define GS_SIZE 4

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
		boost::fusion::at_c<0>(data) = p.template get<0>();

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
	Vcluster & vcl = *global_v_cluster;

	//! Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	if(vcl.getProcessingUnits() != 4)
	return;

	//! Cartesian grid
	size_t sz[2] = { GS_SIZE, GS_SIZE };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph factory
	DistGraphFactory<2, DistGraph_CSR<vx, ed>> g_factory;

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd = g_factory.template construct<NO_EDGE, float, 2 - 1, 0, 1, 2>(sz, box);

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
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 0);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 1);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(2), 2);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(3), 3);
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
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 4);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 5);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(2), 6);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(3), 7);
	}

	if(vcl.getProcessUnitID() == 2)
	{
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 8);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 9);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(2), 10);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(3), 11);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(4), 12);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(5), 13);
	}

	if(vcl.getProcessUnitID() == 3)
	{
		BOOST_REQUIRE_EQUAL(gd.getVertexId(0), 14);
		BOOST_REQUIRE_EQUAL(gd.getVertexId(1), 15);
	}
}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_redistribution)
{
	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

	if(vcl.getProcessingUnits() != 4)
	return;

	//! Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! Cartesian grid
	size_t sz[2] = { 4, 4 };

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 1.0, 1.0 });

	//! Distributed graph factory
	DistGraphFactory<2, DistGraph_CSR<vx, ed>> g_factory;

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd = g_factory.template construct<NO_EDGE, float, 2 - 1, 0, 1, 2>(sz, box);

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

	VTKWriter<DistGraph_CSR<vx, ed>, DIST_GRAPH> gv2(gd);
	gv2.write("dist_graph_redistribution_0.vtk");

	if(vcl.getProcessUnitID() == 0)
	{
		bool test = compare("dist_graph_redistribution_0.vtk", "dist_graph_redistribution_0_test.vtk");
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
		bool test = compare("dist_graph_redistribution_1.vtk","dist_graph_redistribution_1_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}

}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_cartesian)
{
	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

//	if(vcl.getProcessingUnits() != 3)
//	return;

	// non-periodic boundary condition
	size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	//! [Create CartDecomposition]
	CartDecomposition<3, float> dec(vcl);

	// Physical domain
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });
	size_t div[3] = {8,8,8};

	// Grid size and info
	size_t gsz[3] = {8,8,8};
	grid_sm<3,void> g_sm(gsz);

	// Define ghost
	Ghost<3, float> g(0.01);

	// Decompose
	dec.setParameters(div, box, bc, g);
	dec.decompose();

	grid_dist_id_iterator_dec<CartDecomposition<3,float>> it_dec(dec,gsz);

	size_t cnt = 0;
	while (it_dec.isNext())
	{
		cnt++;
		++it_dec;
	}

	openfpm::vector<size_t> v_cnt(vcl.getProcessingUnits());

	// Sent and receive the size of each subgraph
	vcl.allGather(cnt, v_cnt);
	vcl.execute();

	cnt = 0;
	for (long int i = 0; i <= ((long int)vcl.getProcessUnitID()) - 1 ; ++i)
		cnt += v_cnt.get(i);

	// count the points

	//! Distributed graph
	DistGraph_CSR<aggregate<size_t[3]>, aggregate<size_t>> dg;

	grid_dist_id_iterator_dec<CartDecomposition<3,float>> it_dec2(dec,gsz);
	while (it_dec2.isNext())
	{
		auto key = it_dec2.get();

		aggregate<size_t[3]> v;

		v.template get<0>()[0] = key.get(0);
		v.template get<0>()[1] = key.get(1);
		v.template get<0>()[2] = key.get(2);

		size_t gid = g_sm.LinId(key);
		dg.add_vertex(v, gid, cnt);

		cnt++;
		++it_dec2;
	}

	dg.initProperties();

	// we ask for some random vertex

	std::default_random_engine rg;
	std::uniform_int_distribution<size_t> d(0,g_sm.size()-1);

	openfpm::vector<size_t> v_req;

/*	for (size_t i = 0 ; i < 16 ; i++)
	{
		size_t v = d(rg);*/

		if (vcl.getProcessUnitID() == 0)
			dg.reqVertex(450);

/*		dg.reqVertex(v);
	}*/

	dg.sync();

	if (vcl.getProcessUnitID() == 0)
	{
		grid_key_dx<3> key;
		// get the position information
		key.set_d(0,dg.getVertex(450).template get<0>()[0]);
		key.set_d(1,dg.getVertex(450).template get<0>()[1]);
		key.set_d(2,dg.getVertex(450).template get<0>()[2]);

		size_t lin_id = g_sm.LinId(key);

	//	BOOST_REQUIRE_EQUAL(lin_id,v_req.get(i));

		std::cout << "Error: " << "   " << lin_id << "    " << key.to_string() << "\n";
	}

/*	for (size_t i = 0 ; i < 16 ; i++)
	{
		grid_key_dx<3> key;
		// get the position information
		key.set_d(0,dg.getVertex(v_req.get(i)).template get<0>()[0]);
		key.set_d(1,dg.getVertex(v_req.get(i)).template get<0>()[1]);
		key.set_d(2,dg.getVertex(v_req.get(i)).template get<0>()[2]);

		size_t lin_id = g_sm.LinId(key);

	//	BOOST_REQUIRE_EQUAL(lin_id,v_req.get(i));

		std::cout << "Error: " << i << "   " << lin_id << "  " << v_req.get(i) << "\n";
	}*/

/*	if (vcl.getProcessUnitID() == 0)
		std::cout << "Error: " << i << "   " << lin_id << "  " << v_req.get(i) << "\n";*/
}

BOOST_AUTO_TEST_CASE( dist_map_graph_use_free_add)
{
	//! Vcluster
	Vcluster & vcl = *global_v_cluster;

	if(vcl.getProcessingUnits() != 4)
	return;

	//! Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	//! Box
	Box<2, float> box( { 0.0, 0.0 }, { 10.0, 10.0 });

	//! Distributed graph
	DistGraph_CSR<vx, ed> gd;

	for (size_t i = 0; i < 4; ++i)
	{
		vx v;
		v.template get<vx::x>()[0] = vcl.getProcessUnitID();
		v.template get<vx::x>()[1] = i;
		v.template get<vx::x>()[2] = 0;
		size_t id = vcl.getProcessUnitID()*4 + i;
		gd.add_vertex(v, id, id);
	}

	gd.initProperties();

	gd.add_edge(0,1);
	gd.add_edge(1,2);
	gd.add_edge(2,3);
	gd.add_edge(0,4);
	gd.add_edge(1,5);
	gd.add_edge(2,6);
	gd.add_edge(3,7);
	gd.add_edge(4,5);
	gd.add_edge(5,6);
	gd.add_edge(6,7);
	gd.add_edge(4,8);
	gd.add_edge(5,9);
	gd.add_edge(6,10);
	gd.add_edge(7,11);
	gd.add_edge(8,9);
	gd.add_edge(9,10);
	gd.add_edge(10,11);
	gd.add_edge(8,12);
	gd.add_edge(9,13);
	gd.add_edge(10,14);
	gd.add_edge(11,15);
	gd.add_edge(12,13);
	gd.add_edge(13,14);
	gd.add_edge(14,15);

	gd.syncEdge();

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
		bool test = compare("dist_graph_free_0.vtk", "dist_graph_free_0_test.vtk");
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
		bool test = compare("dist_graph_free_1.vtk", "dist_graph_free_1_test.vtk");
		BOOST_REQUIRE_EQUAL(true,test);
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif
