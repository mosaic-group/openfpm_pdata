/*
 * vector_dist_NN_tests.hpp
 *
 *  Created on: Aug 20, 2016
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "VCluster/VCluster.hpp"
#include "Vector/vector_dist.hpp"

extern void print_test_v(std::string test, size_t sz);

template<unsigned int opt> using VERLET_MEMFAST_OPT = VERLET_MEMFAST<3,float,opt>;
template<unsigned int opt> using VERLET_MEMBAL_OPT = VERLET_MEMBAL<3,float,opt>;
template<unsigned int opt> using VERLET_MEMMW_OPT = VERLET_MEMMW<3,float,opt>;

template<template <unsigned int> class VerletList>
void test_full_nn(long int k)
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	long int big_step = k / 4;
	big_step = (big_step == 0)?1:big_step;

	print_test_v("Testing 3D full NN search k=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D full NN search k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	for (double r_cut = 0.1 ; r_cut < 1.0; r_cut += 0.1)
	{
		// ghost
		Ghost<3,float> ghost(r_cut*1.001);

		typedef  aggregate<float> part_prop;

		// Distributed vector
		vector_dist<3,float, part_prop > vd(k,box,bc,ghost);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// Fill some properties randomly

			vd.getProp<0>(key) = 0.0;

			++it;
		}

		vd.map();

		// sync the ghost
		vd.ghost_get<0>();

		openfpm::vector<openfpm::vector<size_t>> list_idx;
		openfpm::vector<openfpm::vector<size_t>> list_idx2;

		list_idx.resize(vd.size_local());
		list_idx2.resize(vd.size_local());

		for (size_t i = 0 ; i < vd.size_local() ; i++)
		{
			Point<3,float> p = vd.getPos(i);

			for (size_t j = 0 ; j < vd.size_local_with_ghost(); j++)
			{
				Point<3,float> q = vd.getPos(j);

				if (p.distance2(q) < r_cut * r_cut)
					list_idx.get(i).add(j);
			}

			list_idx.get(i).sort();
		}

		auto NN = vd.getCellList(r_cut);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NN.getNNIterator(NN.getCell(xp));

			while (Np.isNext())
			{
				auto q = Np.get();

				Point<3,float> xq = vd.getPos(q);

				if (xp.distance2(xq) < r_cut * r_cut)
					list_idx2.get(it.get().getKey()).add(q);

				++Np;
			}

			list_idx2.get(it.get().getKey()).sort();

			++it;
		}

		BOOST_REQUIRE_EQUAL(list_idx.size(),list_idx2.size());

		bool ret;
		for (size_t i = 0 ; i < list_idx.size() ; i++)
		{
			BOOST_REQUIRE_EQUAL(list_idx.get(i).size(),list_idx2.get(i).size());

			ret = true;
			for (size_t j = 0 ; j < list_idx.get(i).size() ; j++)
				ret &= list_idx.get(i).get(j) == list_idx2.get(i).get(j);

			BOOST_REQUIRE_EQUAL(ret,true);
		}

		///////////////////////////////////

		auto NNv = vd.template getVerlet<VerletList<VL_NON_SYMMETRIC>>(r_cut*1.0001);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NNv.getNNIterator(it.get().getKey());

			list_idx2.get(it.get().getKey()).clear();

			while (Np.isNext())
			{
				auto q = Np.get();

				Point<3,float> xq = vd.getPos(q);

				if (xp.distance2(xq) < r_cut * r_cut)
					list_idx2.get(it.get().getKey()).add(q);

				++Np;
			}

			list_idx2.get(it.get().getKey()).sort();

			++it;
		}

		BOOST_REQUIRE_EQUAL(list_idx.size(),list_idx2.size());

		for (size_t i = 0 ; i < list_idx.size() ; i++)
		{
			BOOST_REQUIRE_EQUAL(list_idx.get(i).size(),list_idx2.get(i).size());

			ret = true;
			for (size_t j = 0 ; j < list_idx.get(i).size() ; j++)
				ret &= list_idx.get(i).get(j) == list_idx2.get(i).get(j);

			BOOST_REQUIRE_EQUAL(ret,true);
		}

		//////////////////////////////////////////

		NNv.clear();
		vd.updateVerlet(NNv,r_cut*1.0001);

		it = vd.getDomainIterator();

		while (it.isNext())
		{
			Point<3,float> xp = vd.getPos(it.get());
			auto Np = NNv.getNNIterator(it.get().getKey());

			list_idx2.get(it.get().getKey()).clear();

			while (Np.isNext())
			{
				auto q = Np.get();

				Point<3,float> xq = vd.getPos(q);

				if (xp.distance2(xq) < r_cut * r_cut)
					list_idx2.get(it.get().getKey()).add(q);

				++Np;
			}

			list_idx2.get(it.get().getKey()).sort();

			++it;
		}

		BOOST_REQUIRE_EQUAL(list_idx.size(),list_idx2.size());

		for (size_t i = 0 ; i < list_idx.size() ; i++)
		{
			BOOST_REQUIRE_EQUAL(list_idx.get(i).size(),list_idx2.get(i).size());

			ret = true;
			for (size_t j = 0 ; j < list_idx.get(i).size() ; j++)
				ret &= list_idx.get(i).get(j) == list_idx2.get(i).get(j);

			BOOST_REQUIRE_EQUAL(ret,true);
		}
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_full_NN )
{
	auto & v_cl = create_vcluster();

#ifdef TEST_COVERAGE_MODE
    long int k = 50 * v_cl.getProcessingUnits();
#else
    long int k = 750 * v_cl.getProcessingUnits();
#endif

	test_full_nn<VERLET_MEMFAST_OPT>(k);

	k /= 2;
	test_full_nn<VERLET_MEMBAL_OPT>(k);
	k /= 2;
	test_full_nn<VERLET_MEMMW_OPT>(k);
}

BOOST_AUTO_TEST_CASE( vector_dist_particle_iteration )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 750 * v_cl.getProcessingUnits();

	print_test_v("Testing 3D particle cell iterator=",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D full NN search k=" << k );

	Box<3,float> box({0.0,0.0,0.0},{1.0,1.0,1.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 0.1;

	// ghost
	Ghost<3,float> ghost(r_cut);

	typedef  aggregate<float> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost,BIND_DEC_TO_GHOST);

		auto it = vd.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			vd.getPos(key)[0] = ud(eg);
			vd.getPos(key)[1] = ud(eg);
			vd.getPos(key)[2] = ud(eg);

			// Fill some properties randomly

			vd.getProp<0>(key) = 0.0;

			++it;
		}

	vd.map();

	// sync the ghost
	vd.ghost_get<0>();

	openfpm::vector<size_t> ids;
	ids.resize(vd.size_local());

	auto NN = vd.getCellListSym(r_cut);

	auto it_pcell = vd.getDomainIteratorCells(NN);

	size_t count = 0;
	while (it_pcell.isNext())
	{
		count++;

		size_t id = it_pcell.get();
		ids.get(id) = 1;

		BOOST_REQUIRE(id < vd.size_local());

		++it_pcell;
	}

	v_cl.sum(count);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL((long int)count,k);
}

BOOST_AUTO_TEST_CASE( vector_dist_particle_NN_update_with_limit )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 750 * v_cl.getProcessingUnits();

	print_test_v("Testing 3D particle cell-list with radius at limit= ",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D particle cell-list with radius at limit= " << k );

	Box<3,float> box({0.0,0.0,0.0},{0.1,0.39,0.39});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 0.1;

	// ghost
	Ghost<3,float> ghost(r_cut);

	typedef  aggregate<float> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ud(eg);
		vd.getPos(key)[1] = ud(eg);
		vd.getPos(key)[2] = ud(eg);

		// Fill some properties randomly

		vd.getProp<0>(key) = 0.0;

		++it;
	}

	vd.map();

	// sync the ghost
	vd.ghost_get<0>();

	auto NN = vd.getCellListSym(r_cut);

	auto cell1 = NN.getCellBox();

	vd.getDecomposition().decompose();
	vd.map();

	vd.updateCellList(NN);

	auto cell2 = NN.getCellBox();

	BOOST_REQUIRE_EQUAL(cell1.getHigh(0),cell2.getHigh(0));
	BOOST_REQUIRE_EQUAL(cell1.getHigh(1),cell2.getHigh(1));
	BOOST_REQUIRE_EQUAL(cell1.getHigh(2),cell2.getHigh(2));
}

BOOST_AUTO_TEST_CASE( vector_dist_particle_getCellListSym_with_div )
{
	Vcluster<> & v_cl = create_vcluster();

	if (v_cl.getProcessingUnits() > 12)
		return;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

    long int k = 750 * v_cl.getProcessingUnits();

	print_test_v("Testing 3D particle getCellListSym with div =",k);
	BOOST_TEST_CHECKPOINT( "Testing 3D particle getCellListSym with div = " << k );

	Box<3,float> box({0.0,0.0,0.0},{0.1,0.39,0.39});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 0.1;

	// ghost
	Ghost<3,float> ghost(r_cut);

	typedef  aggregate<float> part_prop;

	// Distributed vector
	vector_dist<3,float, part_prop > vd(k,box,bc,ghost);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.getPos(key)[0] = ud(eg);
		vd.getPos(key)[1] = ud(eg);
		vd.getPos(key)[2] = ud(eg);

		// Fill some properties randomly

		vd.getProp<0>(key) = 0.0;

		++it;
	}

	vd.map();

	// sync the ghost
	vd.ghost_get<0>();

	auto NN1 = vd.getCellListSym(r_cut);

	size_t div_wp[3] = {NN1.getDivWP()[0],NN1.getDivWP()[1],NN1.getDivWP()[2]};
	size_t pad[3] = {NN1.getPadding()[0],NN1.getPadding()[1],NN1.getPadding()[2]};

	auto NN2 = vd.getCellListSym(div_wp,pad);

	// Check that the two Cell list are identical

	// grid info
	size_t div[3] = {NN1.getInternalGrid().getSize()[0],
					 NN1.getInternalGrid().getSize()[1],
					 NN1.getInternalGrid().getSize()[2]};

	grid_sm<3,void> g_info(div);

	bool match = true;

	// Create a grid iterator
	grid_key_dx_iterator<3> g_it(g_info);

	while (g_it.isNext())
	{
		size_t cell = g_info.LinId(g_it.get());
		size_t n_ele1 = NN1.getNelements(cell);
		size_t n_ele2 = NN2.getNelements(cell);

		match &= n_ele1 == n_ele2;

		++g_it;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}
