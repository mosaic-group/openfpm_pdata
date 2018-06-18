/*
 * sgrid_dist_id_unit_tests.cpp
 *
 *  Created on: Nov 18, 2017
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "Grid/grid_dist_id.hpp"
#include "Point_test.hpp"

////////////////////////////////////// THEESE TEST ARE BROKEN TO REMPOVE OR FIX ////


const int x = 0;
const int y = 1;
const int z = 2;

BOOST_AUTO_TEST_SUITE( sgrid_dist_id_test )

BOOST_AUTO_TEST_CASE( sgrid_dist_id_basic_test_2D)
{
	periodicity<2> bc = {NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<2,double> domain({-0.3,-0.3},{1.0,1.0});

	// grid size
	size_t sz[2];
	sz[0] = 1024;
	sz[1] = 1024;

	// Ghost
	Ghost<2,double> g(0.01);

	sgrid_dist_id<2,double,Point_test<float>> sg(sz,domain,g,bc);

	// create a grid iterator

	auto it = sg.getGridIterator();

	while(it.isNext())
	{
		auto gkey = it.get();
		auto key = it.get_dist();


		long int sx = gkey.get(0) - 512;
		long int sy = gkey.get(1) - 512;

		if (sx*sx + sy*sy < 128*128)
		{
			sg.template insert<0>(key) = 1.0;
		}

		++it;
	}

	sg.write("sg_test_write");

	bool match = true;
	auto it2 = sg.getGridIterator();

	while(it2.isNext())
	{
		auto gkey = it2.get();
		auto key = it2.get_dist();

		long int sx = gkey.get(0) - 512;
		long int sy = gkey.get(1) - 512;

		if (sx*sx + sy*sy < 128*128)
		{
			match &= (sg.template get<0>(key) == 1.0);
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);

	sg.ghost_get<0>();

	auto & gr = sg.getGridInfo();

	auto it3 = sg.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();
		auto gkey = it3.getGKey(key);

		sg.template insert<0>(key) = gr.LinId(gkey) * gr.LinId(gkey);

		++it3;
	}

	// now we check the stencil

	bool good = true;
	auto it4 = sg.getDomainIterator();

	while (it4.isNext())
	{
		auto key = it4.get();

		double lap;

		// Here we check that all point of the stencil are inside*/

/*		auto key_xm1 = it4.move(key,x,-1);
		auto gkey_xm1 = it4.getGkey(key,x,-1);
		auto key_xp1 = it4.move(key,x,1);
		auto key_ym1 = it4.move(key,y,-1);
		auto key_yp1 = it4.move(key,y,1);
		auto key_zm1 = it4.move(key,z,-1);
		auto key_zp1 = it4.move(key,z,1);


		size_t sx = key_xm1 - 512;
		size_t sy = gkey.get(1) - 512;*/

		lap = sg.template get<0>(key.move(x,1)) + sg.template get<0>(key.move(x,-1)) +
			  sg.template get<0>(key.move(y,1)) + sg.template get<0>(key.move(y,-1)) +
			  4.0*sg.template get<0>(key);

		good &= (lap == 4.0);

		++it4;
	}

	BOOST_REQUIRE_EQUAL(good,true);
}

BOOST_AUTO_TEST_CASE( sgrid_dist_id_basic_test)
{
	periodicity<3> bc = {NON_PERIODIC, NON_PERIODIC, NON_PERIODIC};

	// Domain
	Box<3,double> domain({-0.3,-0.3,-0.3},{1.0,1.0,1.0});

	// grid size
	size_t sz[3];
	sz[0] = 1024;
	sz[1] = 1024;
	sz[2] = 1024;

	// Ghost
	Ghost<3,double> g(0.01);

	sgrid_dist_id<3,double,Point_test<float>> sg(sz,domain,g,bc);

	// create a grid iterator over a bilion point

	auto it = sg.getGridIterator();

	while(it.isNext())
	{
		auto gkey = it.get();
		auto key = it.get_dist();

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 128*128)
		{
			sg.template insert<0>(key) = 1.0;
		}

		++it;
	}

	bool match = true;
	auto it2 = sg.getGridIterator();

	while(it2.isNext())
	{
		auto gkey = it2.get();
		auto key = it2.get_dist();

		size_t sx = gkey.get(0) - 512;
		size_t sy = gkey.get(1) - 512;
		size_t sz = gkey.get(2) - 512;

		if (sx*sx + sy*sy + sz*sz < 128*128)
		{
			match &= (sg.template get<0>(key) == 1.0);
		}

		++it2;
	}

	sg.ghost_get<0>();

	auto & gr = sg.getGridInfo();

	auto it3 = sg.getDomainIterator();

	while (it3.isNext())
	{
		auto key = it3.get();
		auto gkey = it3.getGKey(key);

		sg.template insert<0>(key) = gr.LinId(gkey) * gr.LinId(gkey);

		++it3;
	}

	// now we check the stencil

	bool good = true;
	auto it4 = sg.getDomainIterator();

	while (it4.isNext())
	{
		auto key = it4.get();

		double lap;

		lap = sg.template get<0>(key.move(x,1)) + sg.template get<0>(key.move(x,-1)) +
			  sg.template get<0>(key.move(y,1)) + sg.template get<0>(key.move(y,-1)) +
			  sg.template get<0>(key.move(z,1)) + sg.template get<0>(key.move(z,-1)) +
			  6.0*sg.template get<0>(key);

		good &= (lap == 6.0);

		++it4;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}


BOOST_AUTO_TEST_SUITE_END()
