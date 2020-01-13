/*
 * grid_dist_amr_dist_unit_tests.cpp
 *
 *  Created on: Sep 21, 2017
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK

#include "config.h"
#include <boost/test/unit_test.hpp>
#include "grid_dist_amr.hpp"

BOOST_AUTO_TEST_SUITE( grid_dist_amr_test )

/*! \brief Coarsest levels of the grid
 *
 * \param domain Simulation domain
 * \param coars_g coarsest grid resolution
 * \param n_lvl number of levels
 *
 */
template<typename grid_amr>
void Test3D_amr_create_levels(grid_amr & amr_g, Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot_c = (coars_g - 1)*(coars_g - 1)*(coars_g - 1);
	size_t correct_result = 0;
	size_t correct_result_cell = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i < n_lvl ; i++)
	{
		correct_result += coars_g*coars_g*coars_g;
		correct_result_cell += tot_c*fact;
		coars_g = 2*(coars_g - 1) + 1;
		fact *= 8;
	}

	amr_g.initLevels(n_lvl,g_sz);


	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridIterator(i);

		while (it.isNext())
		{
			auto key = it.get_dist();
			auto akey = amr_g.getAMRKey(i,key);

			amr_g.template insert<0>(akey) = 3.0;

			count++;

			++it;
		}
	}

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		count++;

		++it;
	}

	Vcluster<> & v_cl = create_vcluster();

	v_cl.sum(count);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(count,correct_result);

	auto itc = amr_g.getDomainIteratorCells();

	size_t count_c = 0;

	while (itc.isNext())
	{
		count_c++;

		++itc;
	}

	v_cl.sum(count_c);
	v_cl.execute();

	auto it_level = amr_g.getDomainIteratorCells(3);

	while (it_level.isNext())
	{
		auto key = it_level.get();

		amr_g.template get<0>(3,key);

		++it_level;
	}

	BOOST_REQUIRE_EQUAL(count_c,correct_result_cell);
}



template<unsigned int dim>
inline bool gr_is_inside(const grid_key_dx<dim> & key, const size_t (& sz)[dim])
{
	for (size_t i = 0 ; i < dim ; i++)
	{
		if (key.get(i) >= (long int)sz[i] || key.get(i) < 0)
		{
			return false;
		}
	}

	return true;
}

template <typename grid>
void Test3D_amr_child_parent_get_no_periodic(grid & amr_g, Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	const int x = 0;
	const int y = 1;
	const int z = 2;

	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot = coars_g*coars_g*coars_g;
	size_t correct_result = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i <  n_lvl ; i++)
	{
		correct_result += tot*fact;
		fact *= 8;
	}

	amr_g.initLevels(n_lvl,g_sz);

	//////// Add something /////

	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridIterator(i);

		while (it.isNext())
		{
			auto key = it.get_dist();
			auto akey = amr_g.getAMRKey(i,key);

			amr_g.template insert<0>(akey) = 3.0;

			count++;

			++it;
		}
	}

	////////////////////////////

	std::string test = amr_g.getSpacing(0).toString();

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();
		auto gkey = it.getGKey();

		amr_g.template insert<0>(key) = gkey.get(0);
		amr_g.template insert<1>(key) = gkey.get(1);
		amr_g.template insert<2>(key) = gkey.get(2);

		++it;
	}

	amr_g.template ghost_get<0,1,2>();

	// now we check that move space work

	auto it2 = amr_g.getDomainIterator();

	bool match = true;
	while (it2.isNext())
	{
		auto key = it2.get();
		auto gkey = it2.getGKey();

		auto key_px = key.moveSpace(x,1);
		auto key_gpx = amr_g.getGKey(key_px);

		auto key_mx = key.moveSpace(x,-1);
		auto key_gmx = amr_g.getGKey(key_mx);

		auto key_py = key.moveSpace(y,1);
		auto key_gpy = amr_g.getGKey(key_py);

		auto key_my = key.moveSpace(y,-1);
		auto key_gmy = amr_g.getGKey(key_my);

		auto key_pz = key.moveSpace(z,1);
		auto key_gpz = amr_g.getGKey(key_pz);

		auto key_mz = key.moveSpace(z,-1);
		auto key_gmz = amr_g.getGKey(key_mz);

		if (gr_is_inside(key_gpx,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.template get<0>(key_px) == gkey.get(0) + 1;
			match &= amr_g.template get<1>(key_px) == gkey.get(1);
			match &= amr_g.template get<2>(key_px) == gkey.get(2);


			if (match == false)
			{
				int debug = 0;
				debug++;
			}
		}

		if (gr_is_inside(key_gmx,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.template get<0>(key_mx) == gkey.get(0) - 1;
			match &= amr_g.template get<1>(key_mx) == gkey.get(1);
			match &= amr_g.template get<2>(key_mx) == gkey.get(2);


			if (match == false)
			{
				int debug = 0;
				debug++;
			}
		}

		if (gr_is_inside(key_gpy,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.template get<0>(key_py) == gkey.get(0);
			match &= amr_g.template get<1>(key_py) == gkey.get(1) + 1;
			match &= amr_g.template get<2>(key_py) == gkey.get(2);


			if (match == false)
			{
				int debug = 0;
				debug++;
			}
		}

		if (gr_is_inside(key_gmy,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.template get<0>(key_my) == gkey.get(0);
			match &= amr_g.template get<1>(key_my) == gkey.get(1) - 1;
			match &= amr_g.template get<2>(key_my) == gkey.get(2);


			if (match == false)
			{
				int debug = 0;
				debug++;
			}
		}

		if (gr_is_inside(key_gpz,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.template get<0>(key_pz) == gkey.get(0);
			match &= amr_g.template get<1>(key_pz) == gkey.get(1);
			match &= amr_g.template get<2>(key_pz) == gkey.get(2) + 1;


			if (match == false)
			{
				int debug = 0;
				debug++;
			}
		}

		if (gr_is_inside(key_gmz,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.template get<0>(key_mz) == gkey.get(0);
			match &= amr_g.template get<1>(key_mz) == gkey.get(1);
			match &= amr_g.template get<2>(key_mz) == gkey.get(2) - 1;


			if (match == false)
			{
				int debug = 0;
				debug++;
			}
		}


		// Test to go to all the levels down

		size_t lvl = it2.getLvl();

		if (lvl < amr_g.getNLvl() - 1)
		{
			auto key_l1 = key;
			amr_g.moveLvlDw(key_l1);
			auto key_gl1 = amr_g.getGKey(key_l1);

			for (size_t s = 0 ; s < 3 ; s++)
			{
				match &= key_gl1.get(s) >> 1 == gkey.get(s);
				match &= amr_g.template get<0>(key_l1) == key_gl1.get(0);
				match &= amr_g.template get<1>(key_l1) == key_gl1.get(1);
				match &= amr_g.template get<2>(key_l1) == key_gl1.get(2);
			}
		}

		if (lvl != 0)
		{
			auto key_l1 = key;
			amr_g.moveLvlUp(key_l1);
			auto key_gl1 = amr_g.getGKey(key_l1);

			for (size_t s = 0 ; s < 3 ; s++)
			{
				match &= gkey.get(s) >> 1 == key_gl1.get(s);

				match &= amr_g.template get<0>(key_l1) == key_gl1.get(0);
				match &= amr_g.template get<1>(key_l1) == key_gl1.get(1);
				match &= amr_g.template get<2>(key_l1) == key_gl1.get(2);
			}
		}

		if (match == false)
		{
			int debug = 0;
			debug++;
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}


template <typename grid>
void Test3D_amr_child_parent_get_periodic(grid & amr_g, Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	const int x = 0;
	const int y = 1;
	const int z = 2;

	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot = coars_g*coars_g*coars_g;
	size_t correct_result = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i <  n_lvl ; i++)
	{
		correct_result += tot*fact;
		fact *= 8;
	}

	amr_g.initLevels(n_lvl,g_sz);

	//////// Add something /////

	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridIterator(i);

		while (it.isNext())
		{
			auto key = it.get_dist();
			auto akey = amr_g.getAMRKey(i,key);

			amr_g.template insert<0>(akey) = 3.0;

			count++;

			++it;
		}
	}

	////////////////////////////

	std::string test = amr_g.getSpacing(0).toString();

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();
		auto gkey = it.getGKey();

		amr_g.template insert<0>(key) = gkey.get(0);
		amr_g.template insert<1>(key) = gkey.get(1);
		amr_g.template insert<2>(key) = gkey.get(2);

		++it;
	}

	amr_g.template ghost_get<0,1,2>();

	// now we check that move space work

	auto it2 = amr_g.getDomainIterator();

	bool match = true;
	while (it2.isNext())
	{
		auto key = it2.get();
		auto gkey = it2.getGKey();

		auto key_px = key.moveSpace(x,1);
		auto key_mx = key.moveSpace(x,-1);
		auto key_py = key.moveSpace(y,1);
		auto key_my = key.moveSpace(y,-1);
		auto key_pz = key.moveSpace(z,1);
		auto key_mz = key.moveSpace(z,-1);

		match &= amr_g.template get<0>(key_px) == openfpm::math::positive_modulo(gkey.get(0) + 1,amr_g.getGridInfoVoid(it2.getLvl()).size(0));
		match &= amr_g.template get<1>(key_px) == gkey.get(1);
		match &= amr_g.template get<2>(key_px) == gkey.get(2);

		match &= amr_g.template get<0>(key_mx) == openfpm::math::positive_modulo(gkey.get(0) - 1,amr_g.getGridInfoVoid(it2.getLvl()).size(0));
		match &= amr_g.template get<1>(key_mx) == gkey.get(1);
		match &= amr_g.template get<2>(key_mx) == gkey.get(2);

		match &= amr_g.template get<0>(key_py) == gkey.get(0);
		match &= amr_g.template get<1>(key_py) == openfpm::math::positive_modulo(gkey.get(1) + 1,amr_g.getGridInfoVoid(it2.getLvl()).size(1));
		match &= amr_g.template get<2>(key_py) == gkey.get(2);

		match &= amr_g.template get<0>(key_my) == gkey.get(0);
		match &= amr_g.template get<1>(key_my) == openfpm::math::positive_modulo(gkey.get(1) - 1,amr_g.getGridInfoVoid(it2.getLvl()).size(1));
		match &= amr_g.template get<2>(key_my) == gkey.get(2);

		match &= amr_g.template get<0>(key_pz) == gkey.get(0);
		match &= amr_g.template get<1>(key_pz) == gkey.get(1);
		match &= amr_g.template get<2>(key_pz) == openfpm::math::positive_modulo(gkey.get(2) + 1,amr_g.getGridInfoVoid(it2.getLvl()).size(2));

		match &= amr_g.template get<0>(key_mz) == gkey.get(0);
		match &= amr_g.template get<1>(key_mz) == gkey.get(1);
		match &= amr_g.template get<2>(key_mz) == openfpm::math::positive_modulo(gkey.get(2) - 1,amr_g.getGridInfoVoid(it2.getLvl()).size(2));

		// Test to go to all the levels down

		size_t lvl = it2.getLvl();

		if (lvl < amr_g.getNLvl() - 1)
		{
			auto key_l1 = key;
			amr_g.moveLvlDw(key_l1);
			auto key_gl1 = amr_g.getGKey(key_l1);

			for (size_t s = 0 ; s < 3 ; s++)
			{
				match &= key_gl1.get(s) >> 1 == gkey.get(s);
				match &= amr_g.template get<0>(key_l1) == key_gl1.get(0);
				match &= amr_g.template get<1>(key_l1) == key_gl1.get(1);
				match &= amr_g.template get<2>(key_l1) == key_gl1.get(2);
			}
		}

		if (lvl != 0)
		{
			auto key_l1 = key;
			amr_g.moveLvlUp(key_l1);
			auto key_gl1 = amr_g.getGKey(key_l1);

			for (size_t s = 0 ; s < 3 ; s++)
			{
				match &= gkey.get(s) >> 1 == key_gl1.get(s);

				match &= amr_g.template get<0>(key_l1) == key_gl1.get(0);
				match &= amr_g.template get<1>(key_l1) == key_gl1.get(1);
				match &= amr_g.template get<2>(key_l1) == key_gl1.get(2);
			}
		}

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

template <typename grid>
void Test3D_amr_ghost_it(grid & amr_g, Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot = coars_g*coars_g*coars_g;
	size_t correct_result = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i <  n_lvl ; i++)
	{
		correct_result += tot*fact;
		fact *= 8;
	}

	amr_g.initLevels(n_lvl,g_sz);

	//////// Add something /////

	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridGhostIterator(i);

		while (it.isNext())
		{
			auto key = it.get_dist();
			auto gkey = it.get();
			auto akey = amr_g.getAMRKey(i,key);

			amr_g.template insert<0>(akey) = gkey.get(0);
			amr_g.template insert<1>(akey) = gkey.get(1);
			amr_g.template insert<2>(akey) = gkey.get(2);

			if (gkey.get(0) == -1 || gkey.get(0) == (long int)amr_g.getGridInfoVoid(i).size(0) ||
				gkey.get(1) == -1 || gkey.get(1) == (long int)amr_g.getGridInfoVoid(i).size(1) ||
				gkey.get(2) == -1 || gkey.get(2) == (long int)amr_g.getGridInfoVoid(i).size(2))
			{count++;}

			++it;
		}

		size_t tot = (amr_g.getGridInfoVoid(i).size(0) + 2)*
				     (amr_g.getGridInfoVoid(i).size(1) + 2)*
				     (amr_g.getGridInfoVoid(i).size(2) + 2) - amr_g.getGridInfoVoid(i).size();

		auto & v_cl = create_vcluster();

		if (v_cl.size() == 1)
		{
			v_cl.sum(count);
			v_cl.execute();

			BOOST_REQUIRE_EQUAL(tot,count);
		}

		bool match = true;
		auto it2 = amr_g.getDomainIterator(i);

		while (it2.isNext())
		{
			auto key = it2.get();

			// move -x

			auto key_m1 = key.move(0,-1);
			auto key_gm1 = it2.getGKey(key_m1);
			match &= amr_g.template get<0>(i,key_m1) == key_gm1.get(0);

			// move +x

			auto key_p1 = key.move(0,1);
			auto key_gp1 = it2.getGKey(key_p1);
			match &= amr_g.template get<0>(i,key_p1) == key_gp1.get(0);

			// move -y

			key_m1 = key.move(1,-1);
			key_gm1 = it2.getGKey(key_m1);
			match &= amr_g.template get<1>(i,key_m1) == key_gm1.get(1);

			// move +y

			key_p1 = key.move(1,1);
			key_gp1 = it2.getGKey(key_p1);
			match &= amr_g.template get<1>(i,key_p1) == key_gp1.get(1);

			// move -z

			key_m1 = key.move(2,-1);
			key_gm1 = it2.getGKey(key_m1);
			match &= amr_g.template get<2>(i,key_m1) == key_gm1.get(2);

			// move +z

			key_p1 = key.move(2,1);
			key_gp1 = it2.getGKey(key_p1);
			match &= amr_g.template get<2>(i,key_p1) == key_gp1.get(2);

			++it2;
		}

		BOOST_REQUIRE_EQUAL(match,true);
	}
}



template <typename grid>
void Test3D_amr_domain_ghost_it(grid & amr_g, Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot = coars_g*coars_g*coars_g;
	size_t correct_result = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i <  n_lvl ; i++)
	{
		correct_result += tot*fact;
		fact *= 8;
	}

	amr_g.initLevels(n_lvl,g_sz);

	size_t total_all_level = 0;

	//////// Add something /////

	for (size_t i = 0 ; i < amr_g.getNLvl() ; i++)
	{
		// Fill the AMR with something

		size_t count = 0;

		auto it = amr_g.getGridGhostIterator(i);

		while (it.isNext())
		{
			auto key = it.get_dist();
			auto gkey = it.get();
			auto akey = amr_g.getAMRKey(i,key);

			amr_g.template insert<0>(akey) = gkey.get(0);
			amr_g.template insert<1>(akey) = gkey.get(1);
			amr_g.template insert<2>(akey) = gkey.get(2);

			count++;

			++it;
		}

		size_t tot = (amr_g.getGridInfoVoid(i).size(0) + 2)*
				     (amr_g.getGridInfoVoid(i).size(1) + 2)*
				     (amr_g.getGridInfoVoid(i).size(2) + 2);

		auto & v_cl = create_vcluster();

		if (v_cl.size() == 1)
		{
			v_cl.sum(count);
			v_cl.execute();

			BOOST_REQUIRE_EQUAL(tot,count);
		}

		size_t amr_cnt = 0;
		bool match = true;
		auto it2 = amr_g.getDomainGhostIterator(i);

		while (it2.isNext())
		{
			auto key = it2.get();
			auto key_g = it2.getGKey(key);
			match &= amr_g.template get<0>(i,key) == key_g.get(0);

			total_all_level++;
			amr_cnt++;

			++it2;
		}

		BOOST_REQUIRE_EQUAL(amr_cnt,count);
		BOOST_REQUIRE_EQUAL(match,true);
	}

	// test the total iterator

	size_t gtot_count = 0;
	auto tot_it = amr_g.getDomainGhostIterator();

	while (tot_it.isNext())
	{
		gtot_count++;

		++tot_it;
	}

	BOOST_REQUIRE_EQUAL(gtot_count,total_all_level);
}

template<typename grid_amr>
void Test3D_ghost_put(grid_amr & g_dist_amr, long int k)
{
	size_t sz[3] = {(size_t)k,(size_t)k,(size_t)k};

	g_dist_amr.initLevels(4,sz);

	// Grid sm
	grid_sm<3,void> info(sz);

	size_t count = 0;

	for (size_t i = 0 ; i < g_dist_amr.getNLvl() ; i++)
	{
		auto dom = g_dist_amr.getGridIterator(i);

		while (dom.isNext())
		{
			auto key = dom.get_dist();

			g_dist_amr.template insert<0>(i,key) = -6.0;

			// Count the points
			count++;

			++dom;
		}
	}

	// Set to zero the full grid

	{
	auto dom = g_dist_amr.getDomainIterator();

	while (dom.isNext())
	{
		auto key = dom.get();

		g_dist_amr.template insert<0>(key.moveSpace(0,1)) += 1.0;
		g_dist_amr.template insert<0>(key.moveSpace(0,-1)) += 1.0;
		g_dist_amr.template insert<0>(key.moveSpace(1,1)) += 1.0;
		g_dist_amr.template insert<0>(key.moveSpace(1,-1)) += 1.0;
		g_dist_amr.template insert<0>(key.moveSpace(2,1)) += 1.0;
		g_dist_amr.template insert<0>(key.moveSpace(2,-1)) += 1.0;

		++dom;
	}
	}

	bool correct = true;

	// Domain + Ghost iterator
	auto dom_gi = g_dist_amr.getDomainIterator();

	while (dom_gi.isNext())
	{
		auto key = dom_gi.get();

		correct &= (g_dist_amr.template get<0>(key) == 0);

		++dom_gi;
	}

	g_dist_amr.template ghost_put<add_,0>();

	if (count != 0)
	{BOOST_REQUIRE_EQUAL(correct, false);}

	// sync the ghosts
	g_dist_amr.template ghost_get<0>();

	correct = true;

	// Domain + Ghost iterator
	auto dom_gi2 = g_dist_amr.getDomainIterator();

	while (dom_gi2.isNext())
	{
		auto key = dom_gi2.get();

		correct &= (g_dist_amr.template get<0>(key) == 0);

		++dom_gi2;
	}

	BOOST_REQUIRE_EQUAL(correct, true);
}

template <typename> struct Debug;

BOOST_AUTO_TEST_CASE( grid_dist_amr_get_child_test_nop )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Ghost<3,long int> g(1);
	grid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g(domain3,g);

	Test3D_amr_child_parent_get_no_periodic(amr_g,domain3,k,4);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_get_child_test_p )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	periodicity<3> bc = {PERIODIC,PERIODIC,PERIODIC};

	Ghost<3,long int> g(1);
	grid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g(domain3,g,bc);

	Test3D_amr_child_parent_get_periodic(amr_g,domain3,k,4);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Ghost<3,float> g(0.05);
	grid_dist_amr<3,float,aggregate<float>> amr_g(domain3,g);

	Test3D_amr_create_levels(amr_g,domain3,k,4);

	sgrid_dist_amr<3,float,aggregate<float>> amr_g2(domain3,g);

	Test3D_amr_create_levels(amr_g2,domain3,k,4);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_ghost_it_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Ghost<3,long int> g(1);
	grid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g(domain3,g);

	Test3D_amr_ghost_it(amr_g,domain3,k,4);

	sgrid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g2(domain3,g);

	Test3D_amr_ghost_it(amr_g2,domain3,k,4);

	for (size_t i = 0 ; i < amr_g2.getNLvl() ; i++)
	{BOOST_REQUIRE(amr_g2.size_inserted(i) != 0ul);}

	amr_g2.clear();

	for (size_t i = 0 ; i < amr_g2.getNLvl() ; i++)
	{BOOST_REQUIRE_EQUAL(amr_g2.size_inserted(i),0ul);}
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_domain_ghost_it_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Ghost<3,long int> g(1);
	grid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g(domain3,g);

	Test3D_amr_domain_ghost_it(amr_g,domain3,k,4);

	sgrid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g2(domain3,g);

	Test3D_amr_domain_ghost_it(amr_g2,domain3,k,4);

	for (size_t i = 0 ; i < amr_g2.getNLvl() ; i++)
	{BOOST_REQUIRE(amr_g2.size_inserted(i) != 0ul);}

	amr_g2.clear();

	for (size_t i = 0 ; i < amr_g2.getNLvl() ; i++)
	{BOOST_REQUIRE_EQUAL(amr_g2.size_inserted(i),0ul);}
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_get_child_test_low_res )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 2;

	Ghost<3,long int> g(1);
	grid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g(domain3,g);

	Test3D_amr_child_parent_get_no_periodic(amr_g,domain3,k,4);

	sgrid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g2(domain3,g);

	Test3D_amr_child_parent_get_no_periodic(amr_g2,domain3,k,4);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_test_background_value )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	Ghost<3,long int> g(1);

	sgrid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g2(domain3,g);

	size_t g_sz[3] = {4,4,4};

	amr_g2.initLevels(4,g_sz);

	aggregate<long int,long int,long int> bck;
	bck.get<0>() = -57;
	bck.get<1>() = -90;
	bck.get<2>() = -123;

	amr_g2.setBackgroundValue(bck);

	// Get a non existent point to check that
	// the background value work

	grid_dist_key_dx<3> key(0,grid_key_dx<3>({0,0,0}));
	long int bck0 = amr_g2.get<0>(2,key);
	BOOST_REQUIRE_EQUAL(bck0,-57);
	long int bck1 = amr_g2.get<1>(2,key);
	BOOST_REQUIRE_EQUAL(bck1,-90);
	long int bck2 = amr_g2.get<2>(2,key);
	BOOST_REQUIRE_EQUAL(bck2,-123);

	// Now we insert that point and we check the subsequent point
	amr_g2.insert<0>(2,key) = 5;

	grid_dist_key_dx<3> key2(0,grid_key_dx<3>({1,0,0}));
	bck0 = amr_g2.get<0>(2,key2);
	BOOST_REQUIRE_EQUAL(bck0,-57);
	bck1 = amr_g2.get<1>(2,key2);
	BOOST_REQUIRE_EQUAL(bck1,-90);
	bck2 = amr_g2.get<2>(2,key2);
	BOOST_REQUIRE_EQUAL(bck2,-123);

	auto & g_dist_lvl2 = amr_g2.getDistGrid(2);
	g_dist_lvl2.get_loc_grid(0).internal_clear_cache();

	bck0 = amr_g2.get<0>(2,key2);
	BOOST_REQUIRE_EQUAL(bck0,-57);
	bck1 = amr_g2.get<1>(2,key2);
	BOOST_REQUIRE_EQUAL(bck1,-90);
	bck2 = amr_g2.get<2>(2,key2);
	BOOST_REQUIRE_EQUAL(bck2,-123);

}

BOOST_AUTO_TEST_CASE( grid_dist_amr_get_domain_ghost_check )
{
	// Test grid periodic

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int k = 13;

	BOOST_TEST_CHECKPOINT( "Testing grid periodic k<=" << k );

	// Ghost
	Ghost<3,long int> g(1);

	// periodicity
	periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
	grid_dist_amr<3, float, aggregate<long int>> g_dist(domain,g,pr);

	Test3D_ghost_put(g_dist,k);

	// Distributed grid with id decomposition
	sgrid_dist_amr<3, float, aggregate<long int>> sg_dist(domain,g,pr);

	Test3D_ghost_put(sg_dist,k);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_ghost_put_create )
{
	// Test grid periodic

	Box<3,float> domain({0.0,0.0,0.0},{1.0,1.0,1.0});

	Vcluster<> & v_cl = create_vcluster();

	if ( v_cl.getProcessingUnits() > 32 )
	{return;}

	long int k = 13;

	BOOST_TEST_CHECKPOINT( "Testing grid periodic k<=" << k );

	// Ghost
	Ghost<3,long int> g(1);

	// periodicity
	periodicity<3> pr = {{PERIODIC,PERIODIC,PERIODIC}};

	// Distributed grid with id decomposition
	grid_dist_amr<3, float, aggregate<long int>> g_dist(domain,g,pr);

	Test3D_ghost_put(g_dist,k);

	// Distributed grid with id decomposition
	sgrid_dist_amr<3, float, aggregate<long int>> sg_dist(domain,g,pr);

	Test3D_ghost_put(sg_dist,k);
}

BOOST_AUTO_TEST_SUITE_END()
