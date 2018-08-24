/*
 * vector_dist_performance_common.hpp
 *
 *  Created on: Dec 25, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_PERFORMANCE_VECTOR_DIST_PERFORMANCE_COMMON_HPP_
#define SRC_VECTOR_PERFORMANCE_VECTOR_DIST_PERFORMANCE_COMMON_HPP_

#include "Vector/vector_dist.hpp"

/*! \brief Calculate and put particles' forces
 *
 * \param NN Cell list
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> void calc_forces(T & NN, V & vd, float r_cut)
{
	auto it_v = vd.getDomainIterator();

	float sum[dim];

	for (size_t i = 0; i < dim; i++)
		sum[i] = 0;

	while (it_v.isNext())
	{
		//key
		vect_dist_key_dx key = it_v.get();

    	// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

		for (size_t i = 0; i < dim; i++)
			sum[i] = 0;

    	// Get the neighborhood of the particle
    	auto cell_it = NN.template getNNIterator<NO_CHECK>(NN.getCell(p));

    	while(cell_it.isNext())
    	{
    		auto nnp = cell_it.get();

    		// p != q
    		if (nnp == key.getKey())
    		{
    			++cell_it;
    			continue;
    		}

    		Point<dim,float> q = vd.getPos(nnp);

    		if (p.distance2(q) < r_cut*r_cut)
    		{
				//Calculate the forces
    			float num[dim];
    			for (size_t i = 0; i < dim; i++)
    				num[i] = vd.getPos(key)[i] - vd.getPos(nnp)[i];

    			float denom = 0;
    			for (size_t i = 0; i < dim; i++)
					denom += num[i] * num[i];

    			float res[dim];
    			for (size_t i = 0; i < dim; i++)
					res[i] = num[i] / denom;

    			for (size_t i = 0; i < dim; i++)
					sum[i] += res[i];
    		}
			//Next particle in a cell
			++cell_it;
		}

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<prp>(key)[i] += sum[i];

		//Next particle in cell list
		++it_v;
	}
}

/*! \brief For each particle of vd calculate the accumulation of the distances of the neighborhood
 *          particles inside vd2
 *
 *
 * \param NN Cell list vd
 * \param NN2 Cell list vd2
 * \param vd Distributed vector
 * \param vd2 Distributed vector 2
 * \param r_cut Cut-off radius
 *
 */
template<unsigned int dim, unsigned int prp, typename T, typename V> void cross_calc(T & NN, T & NN2, V & vd, V & vd2)
{
	auto it_v = vd.getDomainIterator();

	while (it_v.isNext())
	{
		//key
		vect_dist_key_dx key = it_v.get();

    	// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

    	// Get the neighborhood of the particle
    	auto cell_it = NN2.template getNNIterator<NO_CHECK>(NN2.getCell(p));

    	double sum = 0.0;

    	while(cell_it.isNext())
    	{
    		auto nnp = cell_it.get();

    		Point<dim,float> q = vd2.getPos(nnp);

    		sum += norm(p - q);

			//Next particle in a cell
			++cell_it;
		}

		vd.template getProp<prp>(key) = sum;

		//Next particle in cell list
		++it_v;
	}
}

/*! \brief Initialize a distributed vector
 *
 * \param vd Distributed vector
 * \param v_cl Global vcluster
 * \param k_int Number of particles
 */
template<unsigned int dim, typename v_dist> void vd_initialize(v_dist & vd, Vcluster<> & v_cl, size_t k_int)
{
	// The random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	//! [Create a vector of random elements on each processor 2D]

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
			vd.getPos(key)[i] = ud(eg);

		++it;
	}

	vd.map();
}


/*! \brief Initialize 2 distributed vectors with equally positioned particles
 *
 * \param vd, vd2 Distributed vectors
 * \param v_cl Global vcluster
 * \param k_int Number of particles
 */
template<unsigned int dim, typename v_dist> void vd_initialize_double(v_dist & vd,v_dist & vd2, Vcluster<> & v_cl, size_t k_int)
{
	// The random generator engine
	std::default_random_engine eg(v_cl.getProcessUnitID()*4313);
	std::uniform_real_distribution<float> ud(0.0f, 1.0f);

	//! [Create a vector of random elements on each processor 2D]

	auto it = vd.getIterator();

	while (it.isNext())
	{
		auto key = it.get();

		for (size_t i = 0; i < dim; i++)
		{
			vd.getPos(key)[i] = ud(eg);
			vd2.getPos(key)[i] = vd.getPos(key)[i];
		}

		++it;
	}

	vd.map();
	vd2.map();
}



/*! \brief Calculate and put particles' forces
 *
 * \param NN Cell list hilbert
 * \param vd Distributed vector
 * \param r_cut Cut-off radius
 */
template<unsigned int dim, size_t prp = 0, typename T, typename V> void calc_forces_hilb(T & NN, V & vd, float r_cut)
{
	auto it_cl = NN.getIterator();

	float sum[dim];

	for (size_t i = 0; i < dim; i++)
		sum[i] = 0;

	while (it_cl.isNext())
	{
		//key
		auto key = it_cl.get();

    	// Get the position of the particles
		Point<dim,float> p = vd.getPos(key);

		for (size_t i = 0; i < dim; i++)
			sum[i] = 0;

    	// Get the neighborhood of the particle
    	auto cell_it = NN.template getNNIterator<NO_CHECK>(NN.getCell(p));

    	while(cell_it.isNext())
    	{
    		auto nnp = cell_it.get();

    		// p != q
    		if (nnp == key)
    		{
    			++cell_it;
    			continue;
    		}

    		Point<dim,float> q = vd.getPos(nnp);

    		if (p.distance2(q) < r_cut*r_cut)
    		{
				//Calculate the forces
    			float num[dim];
    			for (size_t i = 0; i < dim; i++)
    				num[i] = vd.getPos(key)[i] - vd.getPos(nnp)[i];

    			float denom = 0;
    			for (size_t i = 0; i < dim; i++)
					denom += num[i] * num[i];

    			float res[dim];
    			for (size_t i = 0; i < dim; i++)
					res[i] = num[i] / denom;

    			for (size_t i = 0; i < dim; i++)
					sum[i] += res[i];
    		}
			//Next particle in a cell
			++cell_it;
		}

		//Put the forces
		for (size_t i = 0; i < dim; i++)
			vd.template getProp<prp>(key)[i] += sum[i];

		//Next particle in cell list
		++it_cl;
	}
}

#endif /* SRC_VECTOR_PERFORMANCE_VECTOR_DIST_PERFORMANCE_COMMON_HPP_ */
