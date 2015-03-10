/*
 * vector_dist_iterator.hpp
 *
 *  Created on: Mar 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_DIST_ITERATOR_HPP_
#define VECTOR_DIST_ITERATOR_HPP_

#include "VCluster.hpp"

template<unsigned int dim, typename device_v>
class vector_dist_iterator
{
	//! vector list counter
	size_t v_c;

	//! List of the grids we are going to iterate
	Vcluster_object_array<device_v> & vList;

	//! Actual iterator
	size_t v_it;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * \param gk std::vector of the local grid
	 *
	 */
	vector_dist_iterator(Vcluster_object_array<device_v> & gk)
	:v_c(0),vList(gk),v_it(0)
	{
	}

	// Destructor
	~vector_dist_iterator()
	{
	}

	/*! \brief operator=
	 *
	 * assign
	 *
	 */
	vector_dist_iterator<dim,device_v> & operator=(const vector_dist_iterator<dim,device_v> & vdi)
	{
		v_c = vdi.v_c;
		vList = vdi.vList;
		v_it = vdi.v_it;

		return *this;
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	vector_dist_iterator<dim,device_v> operator++()
	{
		++v_it;

		// check if a_it is at the end

		if (v_it.isNext() == true)
			return *this;
		else
		{
			// switch to the new grid

			v_c++;

			// get the next grid iterator

			if (v_c < vList.size())
				v_it = vList[v_c].getDomainIterator();
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */

	bool isNext()
	{
		// If there are no other grid stop

		if (v_c >= vList.size())
			return false;

		return true;
	}

	/*! \brief Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	size_t get()
	{
		return vect_dist_key_dx<dim>(v_c,v_it.get());
	}
};


#endif /* VECTOR_DIST_ITERATOR_HPP_ */
