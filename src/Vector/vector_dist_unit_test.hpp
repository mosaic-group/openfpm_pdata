/*
 * vector_dist_unit_test.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_DIST_UNIT_TEST_HPP_
#define VECTOR_DIST_UNIT_TEST_HPP_

#include <random>
#include "Vector/vector_dist.hpp"

BOOST_AUTO_TEST_SUITE( vector_dist_test )

BOOST_AUTO_TEST_CASE( vector_dist_ghost )
{
	// Communication object
	Vcluster & v_cl = *global_v_cluster;

	typedef Point_test<float> p;
	typedef Point<2,float> s;

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	size_t g_div[]= {16,16};

	// processor division on y direction
	size_t point_div = g_div[1] / v_cl.getProcessingUnits();

	// Create a grid info
	grid_sm<2,void> g_info(g_div);

	// Calculate the grid spacing
	Point<2,float> spacing = box.getP2();
	spacing = spacing / g_div;

	// middle spacing
	Point<2,float> m_spacing = spacing / 2;

	// create a sub iterator
	grid_key_dx<2> start(point_div * v_cl.getProcessUnitID(),0);
	grid_key_dx<2> stop(point_div * (v_cl.getProcessUnitID() + 1) - 1,g_div[0]);
	auto g_sub = g_info.getSubIterator(start,stop);

	// Vector of particles
	vector_dist<Point<2,float>, Point_test<float>, Box<2,float>, CartDecomposition<2,float> > vd(g_info.size(),box);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key_v = it.get();
		auto key = g_sub.get();

		// set the particle position

		vd.template getPos<s::x>(key_v)[0] = key.get(0) * spacing[0] + m_spacing[0];
		vd.template getPos<s::x>(key_v)[1] = key.get(1) * spacing[1] + m_spacing[1];

		if (vd.template getPos<s::x>(key_v)[0] >= 1.0 || vd.template getPos<s::x>(key_v)[1] >= 1.0)
		{
			int debug = 0;
			debug++;
		}

		++g_sub;
		++it;
	}

	// Debug write the particles
	vd.write("Particles_before_map.csv");

	// redistribute the particles according to the decomposition
	vd.map();

	// Debug write particles
	vd.write("Particles_after_map.csv");

	// Fill the scalar with the particle position
	const auto & ct = vd.getDecomposition();

	it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// fill with the processor ID where these particle live
		vd.template getProp<p::s>(key) = vd.getPos<s::x>(key)[0] + vd.getPos<s::x>(key)[1] * 16;
		vd.template getProp<p::v>(key)[0] = v_cl.getProcessUnitID();
		vd.template getProp<p::v>(key)[1] = v_cl.getProcessUnitID();
		vd.template getProp<p::v>(key)[2] = v_cl.getProcessUnitID();

		++it;
	}

	// set the ghost based on the radius cut off
	Ghost<2,float> g(spacing.get(0));

	vd.setGhost(g);

	//! Output the decomposition
	ct.write(".");

	// do a ghost get
	vd.template ghost_get<p::s,p::v>();

	// Debug write the particles with GHOST
	vd.write("Particles_with_ghost.csv",WITH_GHOST);

	// Get the decomposition
	const auto & dec = vd.getDecomposition();

	// Get the ghost external boxes
	openfpm::vector<size_t> vb(dec.getNGhostBox());

	// Get the ghost iterator
	auto g_it = vd.getGhostIterator();

	// Check if the ghost particles contain the correct information
	while (g_it.isNext())
	{
		auto key = g_it.get();

		float x0 = vd.getPos<s::x>(key)[0];
		float x1 = vd.getPos<s::x>(key)[1] * 16;
		float scalar = vd.template getProp<p::s>(key);

		// Check the received data
		BOOST_REQUIRE_EQUAL(vd.getPos<s::x>(key)[0] + vd.getPos<s::x>(key)[1] * 16,vd.template getProp<p::s>(key));

		bool is_in = false;
		size_t b = 0;

		// check if the received data is in one of the ghost boxes
		for ( ; b < dec.getNGhostBox() ; b++)
		{
			if (dec.getGhostBox(b).isInside(vd.getPos<s::x>(key)) == true)
			{is_in = true; break;}
		}
		BOOST_REQUIRE_EQUAL(is_in,true);

		// Check that the particle come from the correct processor
		BOOST_REQUIRE_EQUAL(vd.getProp<p::v>(key)[0],dec.getGhostBoxProcessor(b));

		// Add
		vb.get(b)++;

		++g_it;
	}

    CellDecomposer_sm<2,float> cd(SpaceBox<2,float>(box),g_div,0);

	for (size_t i = 0 ; i < vb.size() ; i++)
	{
		// Calculate how many particle should be in the box
		size_t n_point = cd.getGridPoints(dec.getGhostBox(i)).getVolume();

		BOOST_REQUIRE_EQUAL(n_point,vb.get(i));
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_iterator_test_use )
{
	typedef Point_test<float> p;
	typedef Point<2,float> s;

	Vcluster & v_cl = *global_v_cluster;

    // set the seed
	// create the random generator engine
	std::srand(v_cl.getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	Box<2,float> box({0.0,0.0},{1.0,1.0});
	vector_dist<Point<2,float>, Point_test<float>, Box<2,float>, CartDecomposition<2,float> > vd(4096,box);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.template getPos<s::x>(key)[0] = ud(eg);
		vd.template getPos<s::x>(key)[1] = ud(eg);

		++it;
	}

	vd.map();

	// Check if we have all the local particles
	size_t cnt = 0;
	auto & ct = vd.getDecomposition();
	it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		// Check if local
		BOOST_REQUIRE_EQUAL(ct.isLocal(vd.template getPos<s::x>(key)),true);

		cnt++;

		++it;
	}

	//
	v_cl.reduce(cnt);
	v_cl.execute();
	BOOST_REQUIRE_EQUAL(cnt,4096);
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* VECTOR_DIST_UNIT_TEST_HPP_ */
