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
//Number of particles
size_t k = 1000;
// Dimensionality
const size_t dim = 3;

BOOST_AUTO_TEST_CASE( vector_dist_hdf5_save_test )
{
	/////////////////

	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Saving Distributed 3D Vector..." << std::endl;

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
	Ghost<dim,float> ghost(0.1);

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k,box,bc,ghost);

	// Initialize a dist vector
	//vd_initialize<dim>(vd, v_cl, k);

	auto it = vd.getDomainIterator();

	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
		{
			vd.getPos(key)[i] = ud(eg);
			//std::cout << "Value: " << vd.getPos(key)[i] << std::endl;
		}

		++it;
	}

	std::cout << "Size_local: " << vd.size_local_with_ghost() << std::endl;

	vd.map();

	std::cout << "Size_local after map: " << vd.size_local_with_ghost() << std::endl;

	//vd.template ghost_get<0>();

	//std::cout << "Size_local after ghost get: " << vd.size_local_with_ghost() << std::endl;

	auto it_2 = vd.getDomainIterator();

	while (it_2.isNext())
	{
		auto key = it_2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<0>(key)[i] = 0.51234;
		++it_2;
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
		std::cout << "Loading Distributed 3D Vector..." << std::endl;

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
	Ghost<dim,float> ghost(0.1);

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,ghost);

	timer t;
	t.start();
	// Save the vector
    vd.load("vector_dist.h5");
	t.stop();

	std::cout << "Loading time: " << t.getwct() << std::endl;


	size_t n_part = vd.size_local();
	openfpm::vector<size_t> tot_n_part;
	v_cl.allGather(n_part,tot_n_part);
	v_cl.execute();

	size_t sum = 0;

	for (size_t i = 0; i < tot_n_part.size(); i++)
		sum += tot_n_part.get(i);

	// Check total number of real particles
	BOOST_REQUIRE_EQUAL(sum,k);

/*
	auto it = vd.getDomainIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
		{
			std::cout << "Pos: " << vd.getPos(key)[i] << std::endl;
		}

		++it;
	}
*/

	//vd.template ghost_get<0>();


	auto it_2 = vd.getDomainIterator();

	while (it_2.isNext())
	{
		auto key = it_2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			//BOOST_CHECK_CLOSE(vd.template getProp<0>(key)[i],0.51234,0.0001);
			std::cout << "Prop: " << vd.template getProp<0>(key)[i] << std::endl;
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

	const size_t Ng = 128;

	// we create a 128x128x128 Grid iterator
	size_t sz[3] = {Ng,Ng,Ng};

	for (size_t i = 0; i < dim; i++)
		bc[i] = NON_PERIODIC;

	// ghost
	Ghost<dim,float> ghost(1.0/(Ng-2));

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(k,box,bc,ghost);

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
	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
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
	Ghost<dim,float> ghost(1.0/(Ng-2));

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
