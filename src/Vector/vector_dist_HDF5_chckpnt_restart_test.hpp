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

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test )
{
	// Input data
	//Number of particles
	size_t k = 100;
	// Dimensionality
	const size_t dim = 3;

	/////////////////

	Vcluster & v_cl = create_vcluster();

	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << "Saving distributed vector" << std::endl;

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

	// ghost
	Ghost<3,float> ghost(0.5);

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k,box,bc,ghost);

	// Initialize a dist vector
	//vd_initialize<dim>(vd, v_cl, k);

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
			vd.getPos(key)[i] = 0.45123;

		++it;
	}

	vd.map();

	vd.template ghost_get<0>();

	auto it_2 = vd.getIterator();

	while (it.isNext())
	{
		auto key = it_2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<0>(key)[i] = 0.51234;
		++it_2;
	}

	// Save the vector
    vd.save("vector_dist.h5");
}

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_load_test )
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << "Loading distributed vector" << std::endl;

	const size_t dim = 3;

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

	// ghost
	Ghost<3,float> ghost(0.5);

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,ghost);

	vd.load("vector_dist.h5");

	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
			BOOST_CHECK_CLOSE(vd.getPos(key)[i],0.45123,0.0001);

		++it;
	}

	vd.template ghost_get<0>();

	auto it_2 = vd.getIterator();

	while (it.isNext())
	{
		auto key = it_2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			BOOST_CHECK_CLOSE(vd.template getProp<0>(key)[i],0.51234,0.0001);
		++it_2;
	}
}

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test_2 )
{
	// Input data
	// Number of particles
	size_t k = 100;

	//Dimensinality of the space
	const size_t dim = 3;

	/////////////////

	Vcluster & v_cl = create_vcluster();

	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << "Saving distributed vector" << std::endl;


	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	const size_t Ng = 128;

	// we create a 128x128x128 Grid iterator
	size_t sz[3] = {Ng,Ng,Ng};

	for (size_t i = 0; i < dim; i++)
		bc[i] = NON_PERIODIC;

	// ghost
	Ghost<3,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,ghost);

	// Initialize a dist vector
	//vd_initialize<dim>(vd, v_cl, k);

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

	BOOST_REQUIRE_EQUAL(it.getSpacing(0),1.0f/(Ng-1));
	BOOST_REQUIRE_EQUAL(it.getSpacing(1),1.0f/(Ng-1));
	BOOST_REQUIRE_EQUAL(it.getSpacing(2),1.0f/(Ng-1));

	vd.map();

	vd.template ghost_get<0>();

	// The random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	// Create a vector of random elements on each processor

	auto it_2 = vd.getIterator();

	while (it.isNext())
	{
		auto key = it_2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<0>(key)[i] = ud(eg);
			//vd.getPos(key)[i]
		++it_2;
	}

	// Save the vector
    vd.save("vector_dist_2.h5");
}

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_load_test_2 )
{
	if (create_vcluster().getProcessUnitID() == 0)
		std::cout << "Loading distributed vector" << std::endl;

	const size_t dim = 3;

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


	const size_t Ng = 128;

	// we create a 128x128x128 Grid iterator
	size_t sz[3] = {Ng,Ng,Ng};

	// ghost
	Ghost<3,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,ghost);

	vd.load("vector_dist_2.h5");

	auto NN = vd.getCellList(0.5);

	auto it = vd.getGridIterator(sz);

	while (it.isNext())
	{
		auto key = it.get();

		++it;
	}

	BOOST_REQUIRE_EQUAL(it.getSpacing(0),1.0f/(Ng-1));
	BOOST_REQUIRE_EQUAL(it.getSpacing(1),1.0f/(Ng-1));
	BOOST_REQUIRE_EQUAL(it.getSpacing(2),1.0f/(Ng-1));
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* SRC_VECTOR_VECTOR_DIST_HDF5_CHCKPNT_RESTART_TEST_HPP_ */
