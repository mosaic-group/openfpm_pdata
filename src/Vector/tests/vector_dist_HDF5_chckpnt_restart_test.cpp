/*
 * vector_dist_HDF5_save.hpp
 *
 *  Created on: Jun 12, 2016
 *      Author: Yaroslav Zaluzhnyi
 */
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Vector/vector_dist.hpp"
#include "Packer_Unpacker/Pack_selector.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"
#include "Vector/performance/vector_dist_performance_util.hpp"
#include "NN/CellList/CellList_util.hpp"

#include "hdf5.h"

BOOST_AUTO_TEST_SUITE( vd_hdf5_chckpnt_rstrt_test )

// Dimensionality
const size_t dim = 3;

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test )
{
	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	const size_t Ng = 32;

	// we create a Grid iterator
	size_t sz[dim] = {Ng,Ng,Ng};

	for (size_t i = 0; i < dim; i++)
	{bc[i] = NON_PERIODIC;}

	// ghost
	Ghost<dim,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]> > vd(0,box,bc,ghost);

	// Put particles

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		vd.add();

		auto key = it.get();

		vd.getLastPos()[0] = key.get(0) * it.getSpacing(0);
		vd.getLastPos()[1] = key.get(1) * it.getSpacing(1);
		vd.getLastPos()[2] = key.get(2) * it.getSpacing(2);

		++it;
	}

	vd.map();

	// Put forces

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<0>(key)[i] = 0.51234 + vd.getPos(key)[0] + vd.getPos(key)[1]+ vd.getPos(key)[2];

		++it2;
	}

	// Save the vector
    vd.save("vector_dist.h5");

    vector_dist<dim,float, aggregate<float[dim]> > vd2(0,box,bc,ghost);

    vd2.load("vector_dist.h5");

    // Check that vd and vd2 match

    auto it3 = vd.getDomainIterator();

    BOOST_REQUIRE_EQUAL(vd.size_local(),vd2.size_local());

    bool check = true;
	while (it3.isNext())
	{
		auto p = it3.get();

		Point<3,float> p1 = vd.getPos(p);
		Point<3,float> p2 = vd2.getPos(p);

		check &= (p1 == p2);

		check &= (vd.template getProp<0>(p)[0] == vd2.template getProp<0>(p)[0]);
		check &= (vd.template getProp<0>(p)[1] == vd2.template getProp<0>(p)[1]);
		check &= (vd.template getProp<0>(p)[2] == vd2.template getProp<0>(p)[2]);

		++it3;
	}

	BOOST_REQUIRE_EQUAL(check,true);
}



BOOST_AUTO_TEST_CASE( vector_dist_hdf5_load_test )
{
#ifndef SE_CLASS3

	Vcluster & v_cl = create_vcluster();

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	for (size_t i = 0; i < dim; i++)
		bc[i] = NON_PERIODIC;


	const size_t Ng = 32;

	// ghost
	Ghost<dim,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]> > vd(0,box,bc,ghost);

	// Load the vector
    vd.load("test_data/vector_dist_24.h5");

	/////////////////// Checking data ///////////////////////

	// Check total number of particles
	size_t n_part = vd.size_local();
	v_cl.sum(n_part);
	v_cl.execute();

	BOOST_REQUIRE_EQUAL(n_part,Ng*Ng*Ng);

	// Check properties

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		// Check the properties
		for (size_t i = 0; i < dim; i++)
			BOOST_REQUIRE_EQUAL(vd.template getProp<0>(key)[i],(float)(0.51234 + vd.getPos(key)[0] + vd.getPos(key)[1]+ vd.getPos(key)[2]));

		++it2;
	}

#endif
}

BOOST_AUTO_TEST_SUITE_END()

