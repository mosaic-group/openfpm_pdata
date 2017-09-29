/*
 * grid_amr_dist.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: i-bird
 */

#ifndef AMR_GRID_AMR_DIST_HPP_
#define AMR_GRID_AMR_DIST_HPP_

#include "Grid/grid_dist_id.hpp"
#include "Amr/grid_dist_amr_key_iterator.hpp"

#define AMR_IMPL_TRIVIAL 1
#define AMR_IMPL_PATCHES 2
#define AMR_IMPL_OPENVDB 3

template<unsigned int dim, typename St, typename T, unsigned int impl=AMR_IMPL_TRIVIAL ,typename Decomposition = CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_amr
{

};

/*! \brief AMR Adaptive Multi Resolution Grid
 *
 * \tparam dim Dimensionality
 * \tparam St type of space
 * \tparam T what each point of the grid store
 * \tparam Decomposition type of decomposition
 *
 */
template<unsigned int dim, typename St, typename T,typename Decomposition,typename Memory, typename device_grid >
class grid_dist_amr<dim,St,T,AMR_IMPL_TRIVIAL,Decomposition,Memory,device_grid>
{
	//! Simulation domain
	Box<dim,St> domain;

	//! Ghost continuos
	Ghost<dim,St> g;

	//! Ghost integer
	Ghost<dim,long int> g_int;

	//! array of grids
	openfpm::vector<grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>> gd_array;

	//! Iterator for each distributed grid
	openfpm::vector<grid_dist_iterator<dim,device_grid,FREE>> git;

public:

	/*! \brief Constructor
	 *
	 * \param domain Simulation domain
	 * \param g ghost extension
	 *
	 */
	grid_dist_amr(const Box<dim,St> & domain)
	:domain(domain)
	{

	}

	/*! \brief Constructor
	 *
	 * \param domain Simulation domain
	 * \param g ghost extension
	 *
	 */
	grid_dist_amr(const Box<dim,St> & domain, const Ghost<dim,long int> & g)
	:domain(domain)
	{
	}

	/*! \brief Initialize the amr grid
	 *
	 * \paran n_lvl maximum number of levels (0 mean no additional levels)
	 * \param g_sz coarsest grid size on each direction
	 *
	 */
	void initLevels(size_t n_lvl,const size_t (& g_sz)[dim])
	{
		size_t g_sz_lvl[dim];

		for (size_t i = 0; i < dim ; i++)
		{g_sz_lvl[i] = g_sz[i];}

		// Add the coarse level
		gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g_sz,domain,g_int));

		for (size_t i = 0; i < n_lvl - 1 ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++){g_sz_lvl[j] = g_sz_lvl[j]*2;}
			gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g_sz_lvl,domain,g_int));
		}
	}

	/*! \brief Get domain iterator
	 *
	 * \return an iterator over all the grid levels
	 *
	 */
	grid_dist_amr_key_iterator<dim,device_grid> getDomainIterator()
	{
		git.clear();

		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			git.add(gd_array.get(i).getDomainIterator());
		}

		return grid_dist_amr_key_iterator<dim,device_grid>(git);
	}
};


#endif /* AMR_GRID_AMR_DIST_HPP_ */
