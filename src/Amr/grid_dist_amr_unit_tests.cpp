/*
 * grid_dist_amr_dist_unit_tests.cpp
 *
 *  Created on: Sep 21, 2017
 *      Author: i-bird
 */


#define BOOST_TEST_DYN_LINK

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
void Test3D_amr_create_levels(Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	Ghost<3,float> g(0.05);

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


	grid_dist_amr<3,float,aggregate<float>> amr_g(domain,g);

	amr_g.initLevels(n_lvl,g_sz);

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	size_t count = 0;

	while (it.isNext())
	{
		count++;

		++it;
	}

	Vcluster & v_cl = create_vcluster();

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

void Test3D_amr_child_parent_get(Box<3,float> & domain, size_t coars_g, size_t n_lvl)
{
	const int x = 0;
	const int y = 1;
	const int z = 2;

	Ghost<3,long int> g(1);

	size_t g_sz[3] = {coars_g,coars_g,coars_g};

	size_t tot = coars_g*coars_g*coars_g;
	size_t correct_result = 0;
	size_t fact = 1;

	for (size_t i = 0 ; i <  n_lvl ; i++)
	{
		correct_result += tot*fact;
		fact *= 8;
	}


	grid_dist_amr<3,float,aggregate<long int,long int,long int>> amr_g(domain,g);

	amr_g.initLevels(n_lvl,g_sz);

	std::string test = amr_g.getSpacing(0).toString();

	// Iterate across all the levels initialized
	auto it = amr_g.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();
		auto gkey = it.getGKey();

		amr_g.get<0>(key) = gkey.get(0);
		amr_g.get<1>(key) = gkey.get(1);
		amr_g.get<2>(key) = gkey.get(2);

		++it;
	}

	amr_g.ghost_get<0,1,2>();
	amr_g.write("amr_write_test");

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
			match &= amr_g.get<0>(key_px) == gkey.get(0) + 1;
			match &= amr_g.get<1>(key_px) == gkey.get(1);
			match &= amr_g.get<2>(key_px) == gkey.get(2);
		}

		if (gr_is_inside(key_gmx,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.get<0>(key_mx) == gkey.get(0) - 1;
			match &= amr_g.get<1>(key_mx) == gkey.get(1);
			match &= amr_g.get<2>(key_mx) == gkey.get(2);
		}

		if (gr_is_inside(key_gpy,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.get<0>(key_py) == gkey.get(0);
			match &= amr_g.get<1>(key_py) == gkey.get(1) + 1;
			match &= amr_g.get<2>(key_py) == gkey.get(2);
		}

		if (gr_is_inside(key_gmy,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.get<0>(key_my) == gkey.get(0);
			match &= amr_g.get<1>(key_my) == gkey.get(1) - 1;
			match &= amr_g.get<2>(key_my) == gkey.get(2);
		}

		if (gr_is_inside(key_gpz,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.get<0>(key_pz) == gkey.get(0);
			match &= amr_g.get<1>(key_pz) == gkey.get(1);
			match &= amr_g.get<2>(key_pz) == gkey.get(2) + 1;
		}

		if (gr_is_inside(key_gmz,amr_g.getGridInfoVoid(it2.getLvl()).getSize()) == true)
		{
			match &= amr_g.get<0>(key_mz) == gkey.get(0);
			match &= amr_g.get<1>(key_mz) == gkey.get(1);
			match &= amr_g.get<2>(key_mz) == gkey.get(2) - 1;
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

		++it2;
	}

	BOOST_REQUIRE_EQUAL(match,true);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_get_child_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_amr_child_parent_get(domain3,k,4);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_test )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 16*16*16*create_vcluster().getProcessingUnits();
	k = std::pow(k, 1/3.);

	Test3D_amr_create_levels(domain3,k,4);
}

BOOST_AUTO_TEST_CASE( grid_dist_amr_get_child_test_low_res )
{
	// Domain
	Box<3,float> domain3({0.0,0.0,0.0},{1.0,1.0,1.0});

	long int k = 2;

	Test3D_amr_child_parent_get(domain3,k,4);
}

BOOST_AUTO_TEST_SUITE_END()
