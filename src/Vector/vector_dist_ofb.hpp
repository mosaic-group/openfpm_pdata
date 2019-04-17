/*
 * vector_dist_ofb.hpp
 *
 *  Created on: Jan 13, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_OFB_HPP_
#define SRC_VECTOR_VECTOR_DIST_OFB_HPP_

/*! \brief Out of bound policy it detect out of bound particles and decide what to do
 *
 */


struct KillParticle
{
	/*! \brief It decide what to do when the particle is out
	 *
	 * \param pp_id particle id
	 * \param p_id processor id
	 *
	 * \return -1
	 *
	 */
	static long int out(size_t pp_id, size_t p_id)
	{
		return -1;
	}
};

//! Out-of-bound policy kill the particle but print a warning
struct KillParticleWithWarning
{
	/*! \brief It decide what to do when the particle is out
	 *
	 * \param pp_id particle id
	 * \param p_id processor id
	 *
	 * \return -1
	 *
	 */
	static long int out(size_t pp_id, size_t p_id)
	{
		std::cerr << "Warning: " << __FILE__ << ":" << __LINE__ << " out of bound particle detected " << std::endl;

		return -1;
	}
};

//! Out-of-bound policy do nothing
struct Nothing
{
	/*! \brief It decide what to do when the particle is out
	 *
	 * \param pp_id particle id
	 * \param p_id processor id
	 *
	 * \return -1
	 *
	 */
	static long int out(size_t pp_id, size_t p_id)
	{
		return -1;
	}
};

//! Out-of-bound policy kill the program
struct Error
{
	/*! \brief It decide what to do when the particle is out
	 *
	 * \param pp_id particle id
	 * \param p_id processor id
	 *
	 * \return -1
	 *
	 */
	static long int out(size_t pp_id, size_t p_id)
	{
		std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " out of bound particle detected " << std::endl;

		exit(-1);

		return -1;
	}
};

#endif /* SRC_VECTOR_VECTOR_DIST_OFB_HPP_ */
