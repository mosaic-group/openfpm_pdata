/*
 * vector_dist_iterator.hpp
 *
 *  Created on: Mar 10, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_DIST_ITERATOR_HPP_
#define VECTOR_DIST_ITERATOR_HPP_

#include "Vector/vector_dist_key.hpp"
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

	vector_dist_iterator & operator++()
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
	inline vect_dist_key_dx get()
	{
		vect_dist_key_dx v;
		v.setKey(v_it);
		return v;
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


//! Iterator that Iterate across particle indexes
class vector_dist_iterator_subset
{
	//! Actual iterator
	size_t v_it;

	//! end point
	size_t stop;

	const openfpm::vector<aggregate<int>> & pid;

	public:

	/*! \brief Constructor of the distributed grid
	 *
	 * \param start start point
	 * \param stop end point
	 *
	 */
	vector_dist_iterator_subset(size_t start, size_t stop, const openfpm::vector<aggregate<int>> & pid)
	:v_it(start),stop(stop),pid(pid)
	{
	}

	// Destructor
	~vector_dist_iterator_subset()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	vector_dist_iterator_subset & operator++()
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
	inline vect_dist_key_dx get()
	{
		vect_dist_key_dx v;
		v.setKey(pid.template get<0>(v_it));
		return v;
	}

	/*! \brief Get the subset key
	 *
	 * \return the subset key
	 *
	 */
	inline vect_dist_key_dx getSubset()
	{
		vect_dist_key_dx v;
		v.setKey(v_it);
		return v;
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
