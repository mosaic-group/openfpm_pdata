/*
 * grid_dist_id_HDF5_chckpnt_restart_test.hpp
 *
 *  Created on: Nov 9, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Grid/grid_dist_id.hpp"

BOOST_AUTO_TEST_SUITE( gd_hdf5_chckpnt_rstrt_test )

BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_save_test )
{

	// Input data
	size_t k = 2400;

	float ghost_part = 0.0;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	// grid size
	size_t sz[2];
	sz[0] = k;
	sz[1] = k;

	// Ghost
	Ghost<2,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<2, float, aggregate<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

	// get the decomposition
	auto & dec = g_dist.getDecomposition();

	// check the consistency of the decomposition
	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	size_t count = 0;

	auto it = g_dist.getDomainIterator();

	while (it.isNext())
	{
		//key
		auto key = it.get();

		auto keyg = g_dist.getGKey(key);

		g_dist.template get<0>(key) = keyg.get(0);

		++it;
		count++;
	}

	openfpm::vector<size_t> count_total;
	v_cl.allGather(count,count_total);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < count_total.size(); i++)
	{sum += count_total.get(i);}

	timer t;
	t.start();
	// Save the grid
    g_dist.save("grid_dist_id.h5" + std::to_string(v_cl.getProcessingUnits()));
	t.stop();
}

BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_load_test )
{

	// Input data
	size_t k = 2400;

	float ghost_part = 0.0;

	// Domain
	Box<2,float> domain({0.0,0.0},{1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() >= 32)
		return;

	// grid size
	size_t sz[2];
	sz[0] = k;
	sz[1] = k;

	// Ghost
	Ghost<2,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<2, float, aggregate<float>, CartDecomposition<2,float>> g_dist(sz,domain,g);

	g_dist.load("grid_dist_id.h5" + std::to_string(v_cl.getProcessingUnits()));

	auto it = g_dist.getDomainIterator();

	size_t count = 0;

	bool match = true;
	while (it.isNext())
	{
		//key
		auto key = it.get();

		//BOOST_CHECK_CLOSE(g_dist.template get<0>(key),1,0.0001);
		//std::cout << "Element: " << g_dist.template get<0>(key) << std::endl;

		auto keyg = g_dist.getGKey(key);

		match &= g_dist.template get<0>(key) == keyg.get(0);

		++it;
		count++;
	}

	openfpm::vector<size_t> count_total;
	v_cl.allGather(count,count_total);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < count_total.size(); i++)
		sum += count_total.get(i);

	BOOST_REQUIRE_EQUAL(sum, (size_t)k*k);
	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_2GB_save_test )
{
	float ghost_part = 0.0;

	// Domain
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() != 1)
		return;

	// grid size
	size_t sz[3];
	sz[0] = 970;
	sz[1] = 650;
	sz[2] = 512;

	// Ghost
	Ghost<3,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<3, float, aggregate<double>> g_dist(sz,domain,g);

	// get the decomposition
	auto & dec = g_dist.getDecomposition();

	// check the consistency of the decomposition
	bool val = dec.check_consistency();
	BOOST_REQUIRE_EQUAL(val,true);

	size_t count = 0;

	auto it = g_dist.getDomainIterator();

	while (it.isNext())
	{
		//key
		auto key = it.get();

		auto keyg = g_dist.getGKey(key);

		g_dist.template get<0>(key) = keyg.get(0);

		++it;
		count++;
	}

	openfpm::vector<size_t> count_total;
	v_cl.allGather(count,count_total);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < count_total.size(); i++)
	{sum += count_total.get(i);}

	timer t;
	t.start();
	// Save the grid
    g_dist.save("grid_dist_2GB_id.h5" + std::to_string(v_cl.getProcessingUnits()));
	t.stop();
}

BOOST_AUTO_TEST_CASE( grid_dist_id_hdf5_2GB_load_test )
{
	float ghost_part = 0.0;

	// Domain
	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1,0});

	Vcluster<> & v_cl = create_vcluster();

	// Skip this test on big scale
	if (v_cl.getProcessingUnits() != 1)
		return;

	// grid size
	size_t sz[3];
	sz[0] = 970;
	sz[1] = 650;
	sz[2] = 512;

	// Ghost
	Ghost<3,float> g(ghost_part);

	// Distributed grid with id decomposition
	grid_dist_id<3, float, aggregate<double>> g_dist(sz,domain,g);

	g_dist.load("grid_dist_2GB_id.h5" + std::to_string(v_cl.getProcessingUnits()));

	auto it = g_dist.getDomainIterator();

	size_t count = 0;

	bool match = true;
	while (it.isNext())
	{
		//key
		auto key = it.get();

		//BOOST_CHECK_CLOSE(g_dist.template get<0>(key),1,0.0001);
		//std::cout << "Element: " << g_dist.template get<0>(key) << std::endl;

		auto keyg = g_dist.getGKey(key);

		match &= g_dist.template get<0>(key) == keyg.get(0);

		if (match == false)
		{
			int debug = 0;
			debug++;
		}

		++it;
		count++;
	}

	openfpm::vector<size_t> count_total;
	v_cl.allGather(count,count_total);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < count_total.size(); i++)
		sum += count_total.get(i);

	BOOST_REQUIRE_EQUAL(sum, (size_t)970*650*512);
	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_SUITE_END()

