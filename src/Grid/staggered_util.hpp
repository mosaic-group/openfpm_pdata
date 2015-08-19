/*
 * staggered_util.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_UTIL_HPP_
#define SRC_GRID_STAGGERED_UTIL_HPP_

#include "variadic_to_vmpl.hpp"

/*! Meta-function to apply to the vector
 *
 */
template<typename arg0, typename T>
struct F
{
	typedef grid_cpu<arg0::value,T> type;
};

/*! \brief Create staggered data vector
 *
 * \param s_ele boost::fusion::vector of elements
 *
 */
template <typename s_ele>
struct create_stag_data
{
	typedef v_transform<,s_ele>
};

#endif /* SRC_GRID_STAGGERED_UTIL_HPP_ */
