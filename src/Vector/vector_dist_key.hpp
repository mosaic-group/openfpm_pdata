/*
 * vector_dist_key.hpp
 *
 *  Created on: Mar 10, 2015
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_KEY_HPP_
#define VECTOR_DIST_KEY_HPP_



/*! \brief Grid key for a distributed grid
 *
 * Grid key for a distributed grid
 *
 */

class vect_dist_key_dx
{
	//! grid list counter

	size_t v_c;

	//! Local grid iterator

	size_t key;

public:

	/*! \brief Get the local grid
	 *
	 * \return the id of the local grid
	 *
	 */
	size_t getSub()
	{
		return v_c;
	}

	/*! \brief Get the key
	 *
	 * \return the local key
	 *
	 */
	size_t getKey()
	{
		return key;
	}

	vect_dist_key_dx(int v_c, size_t key)
	:v_c(v_c),key(key)
	{
	}
};



#endif /* VECTOR_DIST_KEY_HPP_ */
