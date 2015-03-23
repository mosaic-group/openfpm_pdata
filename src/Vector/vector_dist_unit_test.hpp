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

BOOST_AUTO_TEST_CASE( vector_dist_iterator_test_use)
{
	typedef Point<2,float> s;

    // set the seed
	// create the random generator engine
	std::srand(global_v_cluster->getProcessUnitID());
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

	auto & ct = vd.getDecomposition();
	it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		BOOST_REQUIRE_EQUAL(ct.isLocal(vd.template getPos<s::x>(key)),true);

		++it;
	}
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* VECTOR_DIST_UNIT_TEST_HPP_ */
