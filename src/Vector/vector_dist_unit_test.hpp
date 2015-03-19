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
	Box<2,float> box({0.0,0.0},{1.0,1.0});
	vector_dist<space<2,float>, Point_test<float>, Box<2,float>, CartDecomposition<2,float> > vd(4096,box);

	// randomized
	size_t seed = global_v_cluster->getProcessUnitID();
	srand (time(NULL)+seed);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		vd.template getPos<space<2,float>::x>(key)[0] = (rand()/(double)(RAND_MAX));
		vd.template getPos<space<2,float>::x>(key)[1] = (rand()/(double)(RAND_MAX));

		++it;
	}


}

BOOST_AUTO_TEST_SUITE_END()

#endif /* VECTOR_DIST_UNIT_TEST_HPP_ */
