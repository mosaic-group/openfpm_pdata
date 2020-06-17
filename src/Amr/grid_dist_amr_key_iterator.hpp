/*
 * grid_amr_dist_key_iterator.hpp
 *
 *  Created on: Sep 22, 2017
 *      Author: i-bird
 */

#ifndef SRC_AMR_GRID_DIST_AMR_KEY_ITERATOR_HPP_
#define SRC_AMR_GRID_DIST_AMR_KEY_ITERATOR_HPP_

#include "Vector/map_vector.hpp"
#include "Grid/Iterators/grid_dist_id_iterator.hpp"
#include "grid_dist_amr_key.hpp"

template<unsigned int dim, typename device_grid, typename device_sub_it, typename it_type = grid_dist_iterator<dim,device_grid,device_sub_it,FREE>>
class grid_dist_amr_key_iterator
{
	//! Array of grid iterators
	openfpm::vector<it_type> & git;

	//! actual it type
	struct actual_it
	{
		it_type & it;
	};

	//! Actual distributed grid iterator
	it_type * a_it;

	//! iterator pointer
	size_t g_c;



	/*! \brief from g_c increment g_c until you find a valid grid
	 *
	 */
	void selectValidGrid()
	{
		// When the grid has size 0 potentially all the other informations are garbage
		while (g_c < git.size() && git.get(g_c).isNext() == false ) g_c++;

		// get the next grid iterator
		if (g_c < git.size())
		{
			a_it = &git.get(g_c);
		}
	}

public:

	/*! \brief Constructor
	 *
	 * \param git vector of iterator
	 *
	 */
	grid_dist_amr_key_iterator(openfpm::vector<it_type> & git)
	:git(git),g_c(0)
	{
		a_it = &git.get(0);

		selectValidGrid();
	}


	//! Destructor
	~grid_dist_amr_key_iterator()
	{
	}


	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */
	inline grid_dist_amr_key_iterator<dim,device_grid,device_sub_it,it_type> & operator++()
	{
		++(*a_it);

		// check if a_it is at the end

		if (a_it->isNext() == true)
		{return *this;}
		else
		{
			// switch to the new iterator
			g_c++;

			selectValidGrid();
		}

		return *this;
	}

	/*! \brief Is there a next point
	 *
	 * \return true is there is a next point
	 *
	 */
	inline bool isNext()
	{
		return g_c < git.size();
	}

	/*! \brief Return the actual AMR grid iterator point
	 *
	 *
	 */
	inline grid_dist_amr_key<dim> get()
	{
		return grid_dist_amr_key<dim>(g_c,a_it->get());
	}

	/*! \brief Return the actual global grid position in the AMR struct in global
	 *         coordinates
	 *
	 *
	 */
	inline grid_key_dx<dim> getGKey()
	{
		return git.get(g_c).getGKey(a_it->get());
	}

	/*! \brief Return the level at which we are
	 *
	 *
	 */
	inline size_t getLvl() const
	{
		return g_c;
	}
};


#endif /* SRC_AMR_GRID_DIST_AMR_KEY_ITERATOR_HPP_ */
