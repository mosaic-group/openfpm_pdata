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
	//! Local grid iterator

	size_t key;

public:

	/*! \brief set the key
	 *
	 * \param key the local key
	 *
	 */
	__device__ __host__ inline void setKey(size_t key)
	{
		this->key = key;
	}

	/*! \brief Get the key
	 *
	 * \return the local key
	 *
	 */
	__device__ __host__ inline size_t getKey() const
	{
		return key;
	}

	/*! \brief Get the key
	 *
	 * \return the local key
	 *
	 */
	__device__ __host__ inline operator unsigned int() const
	{
		return key;
	}

	/*! \brief Convert the key into a string message
	 *
	 * \return a string message
	 *
	 */
	std::string to_string()
	{
		std::stringstream ts;

		ts << "x[0]=" << key;

		return ts.str();
	}

	//! Default constructor
	__device__ __host__ inline vect_dist_key_dx()
	{
		/* coverity[uninit_member] */
	}

	__device__ __host__ inline vect_dist_key_dx(size_t key)
	:key(key)
	{
	}
};



#endif /* VECTOR_DIST_KEY_HPP_ */
