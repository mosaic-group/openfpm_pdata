/*
 * staggered_grid.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_GRID_HPP_
#define SRC_GRID_STAGGERED_GRID_HPP_

typedef boost::mpl::vector stag_elements;

/*! \brief Implementation of the staggered grid
 *
 * \param dim Dimensionality od the staggered grid
 * \param ele elements object on each dimensional objects, must be a stag_elements
 *
 *
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+

		In the case of a 2D staggered grid we have 3 (in general dim+1 ) elements

		+ = vertex
		# = edge
		* = volume

        ele = stag_ele<scalar<float>,Point_test<float>,scalar<float>>

        It place a scalar on (*) an object Point_test<float> on (#) and an object scalar<float> on (+)

 *
 *
 *
 */
template <unsigned int dim, typename ele>
class staggered_grid
{
private:


	openfpm::vector< grid_cpu<dim> >

public:


};


#endif /* SRC_GRID_STAGGERED_GRID_HPP_ */
