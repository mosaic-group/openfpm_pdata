/*
 * vector_dist_iterator.hpp
 *
 *  Created on: Mar 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_DIST_ITERATOR_HPP_
#define VECTOR_DIST_ITERATOR_HPP_

#include "vector_dist_key.hpp"
#include "VCluster/VCluster.hpp"

//! Iterator that Iterate across particle indexes
class vector_dist_iterator
{
	//! Actual iterator
	size_t v_it;

	//! end point
	size_t stop;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * \param start start point
	 * \param stop end point
	 *
	 */
	vector_dist_iterator(size_t start, size_t stop)
	:v_it(start),stop(stop)
	{
	}

	// Destructor
	~vector_dist_iterator()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	vector_dist_iterator operator++()
	{
		++v_it;

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

		if (v_it >= stop)
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
		return vect_dist_key_dx(v_it);
	}

	/*! \brief Reset the iterator
	 *
	 *
	 */
	void reset()
	{
		v_it = 0;
	}
};


#endif /* VECTOR_DIST_ITERATOR_HPP_ */
