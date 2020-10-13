/*
 * ORB_unit_test.hpp
 *
 *  Created on: Mar 16, 2015
 *      Author: i-bird
 */

#ifndef ORB_UNIT_TEST_HPP_
#define ORB_UNIT_TEST_HPP_

#include "ORB.hpp"
#include <random>

BOOST_AUTO_TEST_SUITE( ORB_test )

#define N_POINTS 1024

BOOST_AUTO_TEST_CASE( ORB_test_use)
{
    // set the seed
	// create the random generator engine
	std::srand(create_vcluster().getProcessUnitID());
    std::default_random_engine eg;
    std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	typedef Point<3,float> p;

	// create a random local vector of particles
	openfpm::vector<Point<3,float>> vp(N_POINTS);

	// fill the particles

	auto vp_it = vp.getIterator();

	while (vp_it.isNext())
	{
		auto key = vp_it.get();

		vp.get<p::x>(key)[0] = ud(eg);
		vp.get<p::x>(key)[1] = ud(eg);
		vp.get<p::x>(key)[2] = ud(eg);

		++vp_it;
	}

	// Orthogonal Recursive Bisection
	Box<3,float> dom({0.0,0.0,0.0},{1.0,1.0,1.0});

	ORB<3,float> orb(dom,16,vp);

	//
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* ORB_UNIT_TEST_HPP_ */
