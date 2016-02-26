#ifndef DISTCARTESIAN_GRAPH_UNIT_TEST_HPP
#define DISTCARTESIAN_GRAPH_UNIT_TEST_HPP

#include "Graph/DistCartesianGraphFactory.hpp"
#include "Graph/map_graph.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"
#include "SubdomainGraphNodes.hpp"

#define DGS_SIZE 4

/*!
 *
 * Test node
 *
 */

struct node
{
	//! The node contain 3 unsigned long integer for communication computation and memory
	typedef boost::fusion::vector<size_t, size_t, size_t, size_t, size_t> type;

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
	//! id property id in boost::fusion::vector
	static const unsigned int id = 3;
	static const unsigned int global_id = 4;
	//! total number of properties boost::fusion::vector
	static const unsigned int max_prop = 5;

	node(){

	}

	template<unsigned int id> inline auto get() -> decltype(boost::fusion::at_c < id > (data))
	{
		return boost::fusion::at_c<id>(data);
	}

	template<unsigned int dim, typename Mem> inline node(const encapc<dim, node, Mem> & p)
	{
		this->operator=(p);
	}

	template<unsigned int dim, typename Mem> inline node & operator=(const encapc<dim, node, Mem> & p)
	{
		boost::fusion::at_c<0>(data) = p.template get<0>();
		boost::fusion::at_c<1>(data) = p.template get<1>();
		boost::fusion::at_c<2>(data) = p.template get<2>();
		boost::fusion::at_c<3>(data) = p.template get<3>();
		boost::fusion::at_c<4>(data) = p.template get<4>();

		return *this;
	}

	static bool noPointers()
	{
		return true;
	}
};


const std::string node::attributes::name[] = { "communication", "computation", "memory", "id", "global_id" };

static void * message_receive(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
{
	openfpm::vector<HeapMemory> *v = static_cast<openfpm::vector<HeapMemory> *>(ptr);

	v->get(i).allocate(msg_i);

	return v->get(i).getPointer();

}

BOOST_AUTO_TEST_SUITE (DistCartesianGraphFactory_test)

BOOST_AUTO_TEST_CASE( DistCartesianGraphFactory_3D_use)
{
	// Boundary conditions, non periodic
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Cartesian grid
	size_t sz[3] = { DGS_SIZE, DGS_SIZE, DGS_SIZE };

	// Box
	Box<3, float> box( { 0.0, 0.0, 0.0 }, { 1.0, 1.0, 1.0 });

	// Graph factory
	CartesianGraphFactory<3,Graph_CSR<Point_test<float>,Point_test<float>>> g_factory;

	// Standard graph
	Graph_CSR<Point_test<float>,Point_test<float>> g = g_factory.template construct<NO_EDGE, node::id, float, 3 - 1, 0, 1, 2>(sz, box, bc);

	// Distribution vector
	openfpm::vector<idx_t> vtxdist(vcl.getProcessingUnits() + 1);

	// Distributed graph factory
	DistCartesianGraphFactory<3, DistGraph_CSR<node, node>> dist_g_factory;

	// Distributed graph, communication and computation are just placeholders, not filled coherently
	DistGraph_CSR<node, node> gd = dist_g_factory.template construct<NO_EDGE, node::id, node::global_id, node::communication, node::computation, float, 3 - 1, 0, 1, 2>(sz, box, vtxdist);

	// Check that sizes of graphs are correct
	size_t t_v = gd.getNVertex();
	vcl.sum(t_v);
	vcl.execute();

	BOOST_REQUIRE_EQUAL(g.getNVertex(), t_v);

	size_t count = 0;
	for(size_t i = 0 ; i< g.getNVertex(); i++){
		if(gd.vertexIsInThisGraph(i))
			count++;
	}

	BOOST_REQUIRE_EQUAL(count, gd.getNVertex());

	for(size_t i = (size_t)vtxdist.get(vcl.getProcessUnitID()), local_i = 0; i < (size_t)vtxdist.get(vcl.getProcessUnitID()+1); i++, local_i++)
	{
		BOOST_REQUIRE_EQUAL(gd.vertex(local_i).template get<node::id>(), g.vertex(i).template get<node::id>());
		BOOST_REQUIRE_EQUAL(gd.getNChilds(local_i), g.getNChilds(i));

		for (size_t s = 0; s < g.getNChilds(i); s++)
		{
			BOOST_REQUIRE_EQUAL(gd.getChild(local_i, s), g.getChild(i, s));
		}
	}


}

BOOST_AUTO_TEST_CASE( DistCartesianGraphFactory_2D_use)
{
	// Boundary conditions, non periodic
	size_t bc[] = {NON_PERIODIC,NON_PERIODIC};

	// Vcluster
	Vcluster & vcl = *global_v_cluster;

	if(vcl.getProcessingUnits() != 2)
		return;

	// Initialize the global VCluster
	init_global_v_cluster(&boost::unit_test::framework::master_test_suite().argc,&boost::unit_test::framework::master_test_suite().argv);

	// Cartesian grid
	size_t sz[2] = { DGS_SIZE, DGS_SIZE };

	// Box
	Box<2, float> box( { 0.0, 0.0}, { 1.0, 1.0} );

	// Graph factory
	CartesianGraphFactory<2,Graph_CSR<node, node>> g_factory;

	// Standard graph
	Graph_CSR<node, node> g = g_factory.template construct<NO_EDGE, node::id, float, 2 - 1, 0, 1, 2>(sz, box, bc);

	// Distribution vector
	openfpm::vector<idx_t> vtxdist(vcl.getProcessingUnits() + 1);

	// Distributed graph factory
	DistCartesianGraphFactory<2, DistGraph_CSR<node, node>> dist_g_factory;

	// Distributed graph
	DistGraph_CSR<node, node> gd = dist_g_factory.template construct<NO_EDGE, node::id, node::global_id, node::communication, node::computation, float, 2 - 1, 0, 1, 2>(sz, box, vtxdist);

	// Check that sizes of graphs are correct
	size_t t_v = gd.getNVertex();
	vcl.sum(t_v);
	vcl.execute();

	BOOST_REQUIRE_EQUAL(g.getNVertex(), t_v);

	size_t count = 0;
	for(size_t i = 0 ; i< g.getNVertex(); i++){
		if(gd.vertexIsInThisGraph(i))
			count++;
	}

	BOOST_REQUIRE_EQUAL(count, gd.getNVertex());

	for(size_t i = (size_t)vtxdist.get(vcl.getProcessUnitID()), local_i = 0; i < (size_t)vtxdist.get(vcl.getProcessUnitID()+1); i++, local_i++)
	{
		BOOST_REQUIRE_EQUAL(gd.vertex(local_i).template get<node::id>(), g.vertex(i).template get<node::id>());
		BOOST_REQUIRE_EQUAL(gd.getNChilds(local_i), g.getNChilds(i));

		for (size_t s = 0; s < g.getNChilds(i); s++)
		{
			BOOST_REQUIRE_EQUAL(gd.getChild(local_i, s), g.getChild(i, s));
		}
	}

	if(vcl.getProcessUnitID() == 0){
		std::vector<size_t> pap_prp;
		Packer<size_t,HeapMemory>::packRequest(pap_prp);
		Packer<node,HeapMemory>::packRequest(pap_prp);
		Packer<node,HeapMemory>::packRequest(pap_prp);
		Packer<node,HeapMemory>::packRequest(pap_prp);

		// Calculate how much preallocated memory we need to pack all the objects
		size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

		// allocate the memory
		HeapMemory pmem;
		pmem.allocate(req);
		ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
		mem.incRef();

		Pack_stat sts;

		// try to pack
		Packer<size_t,HeapMemory>::pack(mem, 3, sts);
		Packer<decltype(gd.vertex(0)),HeapMemory>::pack(mem,gd.vertex(0), sts);
		Packer<decltype(gd.vertex(0)),HeapMemory>::pack(mem,gd.vertex(1), sts);
		Packer<decltype(gd.vertex(0)),HeapMemory>::pack(mem,gd.vertex(2), sts);

		vcl.send(1,0,mem.getPointer(0), pmem.size());
		vcl.execute();
	}

	if(vcl.getProcessUnitID() == 1){

		HeapMemory pmem;
		ExtPreAlloc<HeapMemory> mem(sizeof(node)*3 + sizeof(size_t), pmem);

		Unpack_stat ps;
		node v1,v2,v3;

		vcl.recv(0,0,pmem.getPointer(),sizeof(node)*3 + sizeof(size_t));
		vcl.execute();

		size_t size;
		Unpacker<size_t,HeapMemory>::unpack(mem,size,ps);
		BOOST_REQUIRE_EQUAL(size, 3ul);

		Unpacker<node,HeapMemory>::unpack(mem,v1,ps);
		Unpacker<node,HeapMemory>::unpack(mem,v2,ps);
		Unpacker<node,HeapMemory>::unpack(mem,v3,ps);
		BOOST_REQUIRE_EQUAL(v1.template get<node::id>(), 0ul);
		BOOST_REQUIRE_EQUAL(v2.template get<node::id>(), 1ul);
		BOOST_REQUIRE_EQUAL(v3.template get<node::id>(), 2ul);
	}

	//! [Exchange n vertices packed]
	std::vector<size_t> pap_prp;
	Packer<size_t,HeapMemory>::packRequest(pap_prp);
	Packer<node,HeapMemory>::packRequest(pap_prp);
	Packer<node,HeapMemory>::packRequest(pap_prp);
	Packer<node,HeapMemory>::packRequest(pap_prp);

	// Calculate how much preallocated memory we need to pack all the objects
	size_t req = ExtPreAlloc<HeapMemory>::calculateMem(pap_prp);

	// allocate the memory
	HeapMemory pmem;
	pmem.allocate(req);
	ExtPreAlloc<HeapMemory> & mem = *(new ExtPreAlloc<HeapMemory>(pap_prp,pmem));
	mem.incRef();

	Pack_stat sts;

	// try to pack
	Packer<size_t,HeapMemory>::pack(mem, 3, sts);
	for(int i = 0; i < 3; i++)
		Packer<decltype(gd.vertex(0)),HeapMemory>::pack(mem,gd.vertex(i), sts);

	openfpm::vector<size_t> prc;
	openfpm::vector<size_t> size;
	openfpm::vector<void *> ptr;
	openfpm::vector<HeapMemory> packs(vcl.getProcessingUnits());

	if(vcl.getProcessUnitID() == 0)
			prc.add(1);
	if(vcl.getProcessUnitID() == 1)
			prc.add(0);
	size.add(pmem.size());
	ptr.add(mem.getPointer(0));

	// Exchange informations through processors
	vcl.sendrecvMultipleMessagesNBX(prc.size(), &size.get(0), &prc.get(0), &ptr.get(0), message_receive, &packs,
	NONE);

	//! [Exchange n vertices packed]

	if(vcl.getProcessUnitID() == 1){

		Unpack_stat ps;
		openfpm::vector<node> v;

		ExtPreAlloc<HeapMemory> mem(packs.get(0).size(),packs.get(0));

		size_t size;
		Unpacker<size_t,HeapMemory>::unpack(mem,size,ps);

		BOOST_REQUIRE_EQUAL(size, 3ul);

		v.resize(size);

		for(size_t i = 0; i < size; i++){
			node v_n;
			Unpacker<node,HeapMemory>::unpack(mem,v_n,ps);
			v.set(i, v_n);
		}
		BOOST_REQUIRE_EQUAL(v.get(0).template get<node::id>(), 0ul);
		BOOST_REQUIRE_EQUAL(v.get(1).template get<node::id>(), 1ul);
		BOOST_REQUIRE_EQUAL(v.get(2).template get<node::id>(), 2ul);
	}

	if(vcl.getProcessUnitID() == 0){

		Unpack_stat ps;
		openfpm::vector<node> v;

		ExtPreAlloc<HeapMemory> mem(packs.get(1).size(),packs.get(1));

		size_t size;
		Unpacker<size_t,HeapMemory>::unpack(mem,size,ps);

		BOOST_REQUIRE_EQUAL(size, 3ul);

		v.resize(size);

		for(size_t i = 0; i < size; i++){
			node v_n;
			Unpacker<node,HeapMemory>::unpack(mem,v_n,ps);
			v.set(i, v_n);
		}

		BOOST_REQUIRE_EQUAL(v.get(0).template get<node::id>(), 8ul);
		BOOST_REQUIRE_EQUAL(v.get(1).template get<node::id>(), 9ul);
		BOOST_REQUIRE_EQUAL(v.get(2).template get<node::id>(), 10ul);
	}


}

BOOST_AUTO_TEST_SUITE_END()

#endif
