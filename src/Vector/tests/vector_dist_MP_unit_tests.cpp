/*
 * vector_dist_MP_unit_tests.hpp
 *
 *  Created on: Oct 14, 2016
 *      Author: i-bird
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Vector/vector_dist_multiphase_functions.hpp"
#include "VCluster/VCluster.hpp"
#include "Vector/vector_dist.hpp"

BOOST_AUTO_TEST_SUITE( vector_dist_multiphase_test )

BOOST_AUTO_TEST_CASE( vector_dist_multiphase_cell_list_test )
{
	if (create_vcluster().getProcessingUnits() > 24)
		return;

	size_t sz[3] = {60,60,40};

	// The domain
	Box<3,float> box({-1000.0,-1000.0,-1000.0},{2000.0,2000.0,1000.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 51.0;

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	openfpm::vector< vector_dist<3,float, aggregate<double,double>> > phases;

	// first phase
	phases.add( vector_dist<3,float, aggregate<double,double>>(0,box,bc,ghost) );

	// The other 3 phases
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),0) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),0) );
	phases.add( vector_dist<3,float, aggregate<double,double>>(phases.get(0).getDecomposition(),0) );

	// Fill the phases with particles

	auto g_it = phases.get(0).getGridIterator(sz);

	while (g_it.isNext())
	{
		auto key = g_it.get();

		// Add a particle to all the phases
		phases.get(0).add();
		phases.get(1).add();
		phases.get(2).add();
		phases.get(3).add();

		phases.get(0).getLastPos()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(0).getLastPos()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(0).getLastPos()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		phases.get(1).getLastPos()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(1).getLastPos()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(1).getLastPos()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		phases.get(2).getLastPos()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(2).getLastPos()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(2).getLastPos()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		phases.get(3).getLastPos()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(3).getLastPos()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(3).getLastPos()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		++g_it;
	}

	// Sync all phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).map();
	}

	// randomize a little the particles

	for (size_t p = 0 ; p < phases.size() ; p++)
	{
		openfpm::vector<Point<3,float>> vt;

		for (size_t j = 0 ; j < phases.get(p).size_local() ; j++)
		{
			vt.add(phases.get(p).getPos((j + p*133) % phases.get(p).size_local()));
		}
		phases.get(p).getPosVector().swap(vt);
	}

	// Sync all phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).ghost_get<>();
	}

	// Get the cell list of the phase 0 and 1
	auto CL_phase0 = phases.get(0).getCellList(r_cut);
	auto CL_phase1 = phases.get(1).getCellList(r_cut);

	// This function create a Verlet-list between phases 0 and 1
	auto NN_ver01 = createVerlet(phases.get(0),phases.get(1),CL_phase1,r_cut);

	// Check NNver0_1

	bool ret = true;
	auto it = phases.get(0).getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();
		auto Np = NN_ver01.getNNIterator<NO_CHECK>(p.getKey());

		size_t nn_count = 0;

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Count the number of particles
			nn_count++;

			++Np;
		}

		ret &= nn_count == 7ul;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	// Sync all phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).map();
		phases.get(i).ghost_get<>();
	}

	// NN_ver0_all

	// This function create an "Empty" Multiphase Cell List
	auto CL_all = createCellListM<2>(phases,r_cut);

	// This create a Verlet-list between phase 0 and all the other phases
	auto NNver0_all = createVerletM<2>(0,phases.get(0),phases,CL_all,r_cut);

	it = phases.get(0).getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();
		auto Np = NNver0_all.getNNIterator<NO_CHECK>(p.getKey());

		size_t nn_cout[4] = {0,0,0,0};

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Get from which phase it come from
			auto ph_q = Np.getV();

			nn_cout[ph_q]++;

			++Np;
		}

		ret &= nn_cout[0] == 7;
		ret &= nn_cout[1] == 7;
		ret &= nn_cout[2] == 7;
		ret &= nn_cout[3] == 7;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);
}


BOOST_AUTO_TEST_CASE( vector_dist_multiphase_cell_list_sym_test )
{
	if (create_vcluster().getProcessingUnits() > 24)
		return;

	size_t sz[3] = {60,60,40};

	// The domain
	Box<3,float> box({-1000.0,-1000.0,-1000.0},{2000.0,2000.0,1000.0});

	// Boundary conditions
	size_t bc[3]={PERIODIC,PERIODIC,PERIODIC};

	float r_cut = 51.0;

	// ghost, big enough to contain the interaction radius
	Ghost<3,float> ghost(r_cut);

	openfpm::vector< vector_dist<3,float, aggregate<size_t>> > phases;

	// first phase
	phases.add( vector_dist<3,float, aggregate<size_t>>(0,box,bc,ghost) );

	// The other 3 phases
	phases.add( vector_dist<3,float, aggregate<size_t>>(phases.get(0).getDecomposition(),0) );
	phases.add( vector_dist<3,float, aggregate<size_t>>(phases.get(0).getDecomposition(),0) );
	phases.add( vector_dist<3,float, aggregate<size_t>>(phases.get(0).getDecomposition(),0) );

	// Fill the phases with particles

	auto g_it = phases.get(0).getGridIterator(sz);

	while (g_it.isNext())
	{
		auto key = g_it.get();

		// Add a particle to all the phases
		phases.get(0).add();
		phases.get(1).add();
		phases.get(2).add();
		phases.get(3).add();

		phases.get(0).getLastPosWrite()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(0).getLastPosWrite()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(0).getLastPosWrite()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		phases.get(1).getLastPosWrite()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(1).getLastPosWrite()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(1).getLastPosWrite()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		phases.get(2).getLastPosWrite()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(2).getLastPosWrite()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(2).getLastPosWrite()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		phases.get(3).getLastPosWrite()[0] = key.get(0) * g_it.getSpacing(0) + box.getLow(0);
		phases.get(3).getLastPosWrite()[1] = key.get(1) * g_it.getSpacing(1) + box.getLow(1);
		phases.get(3).getLastPosWrite()[2] = key.get(2) * g_it.getSpacing(2) + box.getLow(2);

		++g_it;
	}

	// Sync all phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).map();
	}

	// randomize a little the particles

	for (size_t p = 0 ; p < phases.size() ; p++)
	{
		openfpm::vector<Point<3,float>> vt;

		for (size_t j = 0 ; j < phases.get(p).size_local() ; j++)
		{
			vt.add(phases.get(p).getPos((j + p*133) % phases.get(p).size_local()));
		}
		phases.get(p).getPosVector().swap(vt);
	}

	// Sync all phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).ghost_get<>();
	}

	// Get the cell list of the phase 0 and 1
	auto CL_phase0 = phases.get(0).getCellListSym(r_cut);
	auto CL_phase1 = phases.get(1).getCellListSym(r_cut);

	// This function create a Verlet-list between phases 0 and 1
	auto NN_ver01 = createVerletSym(phases.get(0),phases.get(1),CL_phase1,r_cut);

	// Check NNver0_1

	bool ret = true;
	auto it = phases.get(0).getDomainIterator();

	while (it.isNext())
	{
		auto p = it.get();
		auto Np = NN_ver01.getNNIterator<NO_CHECK>(p.getKey());

		// For each neighborhood of the particle p
		while (Np.isNext())
		{
			// Neighborhood particle q
			auto q = Np.get();

			phases.get(0).getPropWrite<0>(p)++;
			phases.get(1).getPropWrite<0>(q)++;

			++Np;
		}

		++it;
	}

	phases.get(0).ghost_put<add_,0>();
	phases.get(1).ghost_put<add_,0>();

#ifdef SE_CLASS3

	phases.get(1).getDomainIterator();

#endif

	it = phases.get(0).getDomainIterator();
	while (it.isNext())
	{
		auto p = it.get();

		ret &= phases.get(0).getPropRead<0>(p) == 7;
		ret &= phases.get(1).getPropRead<0>(p) == 7;

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);

	// Sync all phases
	for (size_t i = 0 ; i < 4 ; i++)
	{
		phases.get(i).map();
		phases.get(i).ghost_get<>();
	}

	// Reset counter on all phases

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		it = phases.get(i).getDomainAndGhostIterator();
		while (it.isNext())
		{
			auto p = it.get();

			phases.get(i).getPropWrite<0>(p) = 0;

			++it;
		}
	}

	// NN_ver0_all

	// This function create an "Empty" Multiphase Cell List
	auto CL_all = createCellListSymM<2>(phases,r_cut);

	typedef decltype(createVerletSymM<2>(0,phases.get(0),phases,CL_all,r_cut)) verlet_type;

	verlet_type NNver_all[4];

	// This create a Verlet-list between each phase to all the other phases
	NNver_all[0] = createVerletSymM<2>(0,phases.get(0),phases,CL_all,r_cut);
	NNver_all[1] = createVerletSymM<2>(1,phases.get(1),phases,CL_all,r_cut);
	NNver_all[2] = createVerletSymM<2>(2,phases.get(2),phases,CL_all,r_cut);
	NNver_all[3] = createVerletSymM<2>(3,phases.get(3),phases,CL_all,r_cut);

	// all phases to all phases

	for (size_t i = 0 ; i < phases.size() ; i++)
	{
		it = phases.get(i).getDomainIterator();

		while (it.isNext())
		{
			auto p = it.get();
			auto Np = NNver_all[i].getNNIterator<NO_CHECK>(p.getKey());

			// For each neighborhood of the particle p
			while (Np.isNext())
			{
				// Get the particle q near to p
				auto q = Np.getP();

				// Get from which phase it come from
				auto ph_q = Np.getV();

				phases.get(i).getPropWrite<0>(p)++;
				phases.get(ph_q).getPropWrite<0>(q)++;

				++Np;
			}

			++it;
		}
	}

	phases.get(0).ghost_put<add_,0>();
	phases.get(1).ghost_put<add_,0>();
	phases.get(2).ghost_put<add_,0>();
	phases.get(3).ghost_put<add_,0>();

#ifdef SE_CLASS3

	it = phases.get(1).getDomainIterator();
	it = phases.get(2).getDomainIterator();
	it = phases.get(3).getDomainIterator();

#endif

	it = phases.get(0).getDomainIterator();
	while (it.isNext())
	{
		auto p = it.get();

		ret &= phases.get(0).getPropRead<0>(p) == 32;
		ret &= phases.get(1).getPropRead<0>(p) == 32;
		ret &= phases.get(2).getPropRead<0>(p) == 32;
		ret &= phases.get(3).getPropRead<0>(p) == 32;

		if (ret == false)
		{
			std::cout << phases.get(0).getPropRead<0>(p) << std::endl;
			std::cout << phases.get(1).getPropRead<0>(p) << std::endl;
			std::cout << phases.get(2).getPropRead<0>(p) << std::endl;
			std::cout << phases.get(3).getPropRead<0>(p) << std::endl;
		}

		++it;
	}

	BOOST_REQUIRE_EQUAL(ret,true);
}

BOOST_AUTO_TEST_SUITE_END()

