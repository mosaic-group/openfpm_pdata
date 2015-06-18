/*
 * vector_dist_iterator.hpp
 *
 *  Created on: Mar 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_DIST_ITERATOR_HPP_
#define VECTOR_DIST_ITERATOR_HPP_

#include "vector_dist_key.hpp"
#include "VCluster.hpp"

template<typename device_v>
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
	 * \param gk the set of local vectors
	 * \param offset iterator starting point
	 *
	 */
	vector_dist_iterator(Vcluster_object_array<device_v> & gk, size_t offset = 0)
	:v_c(0),vList(gk),v_it(offset)
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
	vector_dist_iterator<device_v> & operator=(const vector_dist_iterator<device_v> & vdi)
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

	vector_dist_iterator<device_v> operator++()
	{
		++v_it;

		// check if a_it is at the end

		if (v_it < vList.get(v_c).size())
			return *this;
		else
		{
			// switch to the new grid

			v_c++;

			// get the next grid iterator

			if (v_c < vList.size())
				v_it = 0;
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
	vect_dist_key_dx get()
	{
		return vect_dist_key_dx(v_c,v_it);
	}
};


#endif /* VECTOR_DIST_ITERATOR_HPP_ */
