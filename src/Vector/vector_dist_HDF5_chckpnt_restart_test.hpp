/*
 * vector_dist_HDF5_save.hpp
 *
 *  Created on: Jun 12, 2016
 *      Author: Yaroslav Zaluzhnyi
 */

#ifndef SRC_VECTOR_VECTOR_DIST_HDF5_CHCKPNT_RESTART_TEST_HPP_
#define SRC_VECTOR_VECTOR_DIST_HDF5_CHCKPNT_RESTART_TEST_HPP_

#include "vector_dist.hpp"
#include "Packer_Unpacker/Pack_selector.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"
#include "vector_dist_performance_util.hpp"


#include "hdf5.h"

BOOST_AUTO_TEST_SUITE( vd_hdf5_chckpnt_rstrt_test )

// Input data

// Number of particles
size_t k = 1000000;

// Dimensionality
const size_t dim = 3;

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Saving distributed vector" << std::endl;

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	const size_t Ng = cbrt(k);

	// we create a Grid iterator
	size_t sz[dim] = {Ng,Ng,Ng};

	for (size_t i = 0; i < dim; i++)
		bc[i] = NON_PERIODIC;

	// ghost
	Ghost<dim,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,ghost);

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

	//BOOST_REQUIRE_EQUAL(it.getSpacing(0),1.0f/(Ng-1));
	//BOOST_REQUIRE_EQUAL(it.getSpacing(1),1.0f/(Ng-1));
	//BOOST_REQUIRE_EQUAL(it.getSpacing(2),1.0f/(Ng-1));

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

	timer t;
	t.start();
	// Save the vector
    vd.save("vector_dist.h5");
	t.stop();

	std::cout << "Saving time: " << t.getwct() << std::endl;
}



BOOST_AUTO_TEST_CASE( vector_dist_hdf5_load_test )
{
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Loading distributed vector" << std::endl;

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


	const size_t Ng = cbrt(k);

	// we create a Grid iterator
	size_t sz[3] = {Ng,Ng,Ng};

	// ghost
	Ghost<dim,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,ghost);

	vd.load("vector_dist.h5");

	timer t;
	t.start();
	// Save the vector
    vd.load("vector_dist.h5");
	t.stop();

	std::cout << "Loading time: " << t.getwct() << std::endl;

	/////////////////// Checking data ///////////////////////

	// Check total number of particles
	size_t n_part = vd.size_local();
	openfpm::vector<size_t> tot_n_part;
	v_cl.allGather(n_part,tot_n_part);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < tot_n_part.size(); i++)
		sum += tot_n_part.get(i);

	BOOST_REQUIRE_EQUAL(sum,k);

	//std::cout << "Sum: " << sum << std::endl;

    // Check spacing (positions)

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		//auto key = it.get();

		++it;
	}

	BOOST_REQUIRE_EQUAL(it.getSpacing(0),1.0f/(Ng-1));
	BOOST_REQUIRE_EQUAL(it.getSpacing(1),1.0f/(Ng-1));
	BOOST_REQUIRE_EQUAL(it.getSpacing(2),1.0f/(Ng-1));


	// Check properties

	auto it2 = vd.getDomainIterator();

	while (it2.isNext())
	{
		auto key = it2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			BOOST_CHECK_CLOSE(vd.template getProp<0>(key)[i],0.51234 + vd.getPos(key)[0] + vd.getPos(key)[1]+ vd.getPos(key)[2],0.0001);

		++it2;
	}
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* SRC_VECTOR_VECTOR_DIST_HDF5_CHCKPNT_RESTART_TEST_HPP_ */
