/*
 * ORB_unit_test.hpp
 *
 *  Created on: Mar 16, 2015
 *      Author: i-bird
 */

#ifndef ORB_UNIT_TEST_HPP_
#define ORB_UNIT_TEST_HPP_

#include "ORB.hpp"
#include "ORBDecompositionStrategy.hpp"
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

/*	openfpm::vector<::Box<3, size_t>> garbage;
	openfpm::vector<SpaceBox<3, float>> sub_domains;
	Box<3,float> bbox;

	orb.convertToSubDomains(garbage,sub_domains,bbox);

	// Create a writer and write
	VTKWriter<decltype(sub_domains),VECTOR_BOX> vtk_box;
	vtk_box.add(sub_domains);
	vtk_box.write("vtk_box.vtk");*/
}

BOOST_AUTO_TEST_CASE( ORB_decomposition_strategy_test_use)
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

	OrbDecompositionStrategy<3,float> orb(create_vcluster());

	const size_t div_[3] = {0,0,0};
	Ghost<3,float> g(0.0);
	const size_t bc[3] = {NON_PERIODIC,NON_PERIODIC,NON_PERIODIC};

	orb.setParameters(div_,dom,bc,g);
	orb.decompose(vp);

	openfpm::vector<::Box<3, size_t>> garbage;
	openfpm::vector<SpaceBox<3, float>> sub_domains;
	Box<3,float> bbox;

	orb.convertToSubDomains(garbage,sub_domains,bbox);

	// Create a writer and write
	VTKWriter<decltype(sub_domains),VECTOR_BOX> vtk_box;
	vtk_box.add(sub_domains);
	vtk_box.write("vtk_box_" + std::to_string(create_vcluster().rank()) + ".vtk");
}

BOOST_AUTO_TEST_SUITE_END()

#endif /* ORB_UNIT_TEST_HPP_ */
