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

	const size_t dim = 3;

	size_t k = 100;

	size_t ghost_part = 0.1;

	/////////////////

	Vcluster & v_cl = create_vcluster();

	if (v_cl.getProcessUnitID() == 0)
		std::cout << "Testing " << dim << "D vector HDF5 save/load" << std::endl;

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	for (size_t i = 0; i < dim; i++)
		bc[i] = PERIODIC;

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,Ghost<dim,float>(ghost_part));

	// Initialize a dist vector
	//vd_initialize<dim>(vd, v_cl, k);

	size_t sz[3] = {10,10,10};

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

	vd.template ghost_get<0>();

	// The random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	//! [Create a vector of random elements on each processor 2D]

	auto it_2 = vd.getIterator();

	while (it.isNext())
	{
		auto key = it_2.get();

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<0>(key)[i] = ud(eg);

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

	size_t ghost_part = 0.1;

	Box<dim,float> box;

	for (size_t i = 0; i < dim; i++)
	{
		box.setLow(i,0.0);
		box.setHigh(i,1.0);
	}

	// Boundary conditions
	size_t bc[dim];

	for (size_t i = 0; i < dim; i++)
		bc[i] = PERIODIC;

	vector_dist<dim,float, aggregate<float[dim]>, CartDecomposition<dim,float> > vd(0,box,bc,Ghost<dim,float>(ghost_part));

	vd.load("vector_dist.h5");
	/*
	auto NN = vd.getCellList(0.5);

	auto it_v = vd.getDomainIterator();

	while (it_v.isNext())
	{
		//key
		vect_dist_key_dx key = it_v.get();

		size_t count = 0;

		// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

		// Get the neighborhood of the particle
		auto cell_it = NN.template getNNIterator<NO_CHECK>(NN.getCell(p));

		while(cell_it.isNext())
		{
			//Next particle in a cell
			++cell_it;
			count++;
		}

		std::cout << "Count: " << count << std::endl;

		//Next particle in cell list
		++it_v;
	}
*/
}

BOOST_AUTO_TEST_SUITE_END()


#endif /* SRC_VECTOR_VECTOR_DIST_HDF5_CHCKPNT_RESTART_TEST_HPP_ */
