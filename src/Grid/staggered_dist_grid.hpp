/*
 * staggered_grid.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_DIST_GRID_HPP_
#define SRC_GRID_STAGGERED_DIST_GRID_HPP_

#include "Grid/grid_dist_id.hpp"
#include "staggered_dist_grid_util.hpp"
#include "VTKWriter.hpp"


/*! \brief Implementation of the staggered grid
 *
 * \param dim Dimensionality of the staggered grid
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
template<unsigned int dim, typename St, typename T, typename Decomposition,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T>>
class staggered_grid_dist : public grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>
{

public:

	staggered_grid_dist(const size_t (& g_sz)[dim], const Box<dim,St> & domain, const Ghost<dim,St> & ghost)
	:grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g_sz,domain,ghost)
	{}

	openfpm::vector<comb<dim>> c_prp[T::max_prop];

	/*! \brief Set the staggered positions
	 *
	 *
	 */
	template<unsigned int p> void setStagPosition(openfpm::vector<comb<dim>> & cmb)
	{
#ifdef SE_CLASS1
		if (mul_extends< boost::mpl::at<ele::type>::type >::ext() != cmb.size())
			std::cerr << __FILE__ << ":" << __LINE << " error properties has " << mul_extends< boost::mpl::at<ele::type>::type >::ext() << " components, but " << cmb.size() << "has been defined \n";
#endif
		c_prp.get(p) = cmb;
	}

	/*! \brief It set all the properties defined to be staggered on the default location
	 *
	 */
	void setDefaultStagPosition()
	{
		// for each properties

		stag_set_position<dim,typename T::type> ssp(c_prp);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(ssp);
	}

	/*! \brief Write a vtk file with the information of the staggered grid
	 *
	 * \param str vtk output file
	 *
	 */
	void write(std::string str)
	{
		// spacing
		Point<dim,St> spacing = grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>::getSpacing();
		spacing = spacing / 2;
		grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>::write(str,c_prp,spacing);
	}
};

#endif /* SRC_GRID_STAGGERED_DIST_GRID_HPP_ */
