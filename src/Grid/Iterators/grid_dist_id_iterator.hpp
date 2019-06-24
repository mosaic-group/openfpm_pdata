/*
 * grid_dist_id_iterator_sub.hpp
 *
 *  Created on: Feb 4, 2015
 *      Author: Pietro Incardona
 */

#ifndef GRID_DIST_ID_ITERATOR_HPP_
#define GRID_DIST_ID_ITERATOR_HPP_

#define FREE 1
#define FIXED 2
#define ITERATION_ISOLATION 4

#include "Grid/grid_dist_key.hpp"
#include "VCluster/VCluster.hpp"
#include "util/GBoxes.hpp"

template<bool is_free>
struct selvg
{
	template<typename a_it_type, typename gdb_ext_type, typename gList_type>
	static inline void call(a_it_type & a_it, gdb_ext_type & gdb_ext, gList_type & gList, size_t & g_c)
	{
		if (gdb_ext.get(g_c).Dbox.isValid() == false)
		{g_c++;}
		else
		{
			a_it.reinitialize(gList.get(g_c).getIterator(gdb_ext.get(g_c).Dbox.getKP1(),gdb_ext.get(g_c).Dbox.getKP2()));
			if (a_it.isNext() == false)	{g_c++;}
		}
	}
};

template<>
struct selvg<false>
{
	template<typename a_it_type, typename gdb_ext_type, typename gList_type>
	static inline void call(a_it_type & a_it, gdb_ext_type & gdb_ext, gList_type & gList, size_t & g_c)
	{
		// Full iterator (no subset)
		a_it.reinitialize(gList.get(g_c).getIterator());
		if (a_it.isNext() == false)	{g_c++;}
	}
};

/*! \brief Distributed grid iterator
 *
 * Iterator across the local elements of the distributed grid
 *
 * \tparam dim dimensionality of the grid
 * \tparam device_grid type of basic grid
 * \tparam stencil it inject the code to calculate stencil offset
 * \tparam sub_iterator it indicate the sub-iterator type of the device_grid
 *
 */
template<unsigned int dim, typename device_grid, typename device_sub_it, int impl, typename stencil  = no_stencil >
class grid_dist_iterator
{
	//! grid list counter
	size_t g_c;

	//! List of the grids we are going to iterate
	const openfpm::vector<device_grid> & gList;

	//! Extension of each grid: domain and ghost + domain
	const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext;

	//! Actual iterator
	device_sub_it a_it;

	//! stop point (is the grid size)
	grid_key_dx<dim> stop;

	/*! \brief from g_c increment g_c until you find a valid grid
	 *
	 */
	void selectValidGrid()
	{
		do
		{
			if (impl == FREE)
			{
				// When the grid has size 0 potentially all the other informations are garbage
				while (g_c < gList.size() && (gList.get(g_c).size() == 0 || gdb_ext.get(g_c).Dbox.isValid() == false ) ) g_c++;
			}
			else
			{
				// When the grid has size 0 potentially all the other informations are garbage
				while (g_c < gList.size() && (gList.get(g_c).size() == 0 || gdb_ext.get(g_c).GDbox.isValid() == false) ) g_c++;
			}

			// get the next grid iterator
			if (g_c < gList.size())
			{
				// Sub iterator are used
/*				if (impl == FREE)
				{
					if (gdb_ext.get(g_c).Dbox.isValid() == false)
					{g_c++;}
					else
					{
						a_it.reinitialize(gList.get(g_c).getIterator(gdb_ext.get(g_c).Dbox.getKP1(),gdb_ext.get(g_c).Dbox.getKP2()));
						if (a_it.isNext() == false)	{g_c++;}
					}
				}
				else
				{
					// Full iterator (no subset)
					a_it.reinitialize(gList.get(g_c).getIterator());
					if (a_it.isNext() == false)	{g_c++;}
				}*/

				selvg<impl == FREE>::call(a_it,gdb_ext,gList,g_c);
			}
		} while (g_c < gList.size() && a_it.isNext() == false);

	}

	public:

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param gk std::vector of the local grid
	 * \param gdb_ext set of local subdomains
	 * \param stop end point
	 *
	 */
	grid_dist_iterator(const openfpm::vector<device_grid> & gk,
					   const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
					   const grid_key_dx<dim> & stop)
	:g_c(0),gList(gk),gdb_ext(gdb_ext),stop(stop)
	{
		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}


	/*! \brief Constructor of the distributed grid iterator with
	 *         stencil support
	 *
	 * \param gk std::vector of the local grid
	 * \param gdb_ext set of local subdomains
	 * \param stop end point
	 * \param stencil_pnt stencil points
	 *
	 */
	grid_dist_iterator(openfpm::vector<device_grid> & gk,
			           const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
					   const grid_key_dx<dim> & stop,
					   const grid_key_dx<dim> (& stencil_pnt)[stencil::nsp])
	:g_c(0),gList(gk),gdb_ext(gdb_ext),a_it(stencil_pnt),stop(stop)
	{
		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}

	//! Copy constructor
	grid_dist_iterator(const grid_dist_iterator<dim,device_grid,device_sub_it,impl,stencil> & g)
	:g_c(g.g_c),gList(g.gList),gdb_ext(g.gdb_ext),a_it(g.a_it),stop(g.stop)
	{}

	//! Copy constructor
	grid_dist_iterator(grid_dist_iterator<dim,device_grid,device_sub_it,impl,stencil> && g)
	:g_c(g.g_c),gList(g.gList),gdb_ext(g.gdb_ext),a_it(g.a_it),stop(g.stop)
	{}

	//! Destructor
	~grid_dist_iterator()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */
	inline grid_dist_iterator<dim,device_grid,device_sub_it,impl,stencil> & operator++()
	{
		++a_it;

		// check if a_it is at the end

		if (a_it.isNext() == true)
			return *this;
		else
		{
			// switch to the new grid
			g_c++;

			selectValidGrid();
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */
	inline bool isNext() const
	{
		// If there are no other grid stop

		if (g_c >= gList.size())
		{return false;}

		return true;
	}

	/*! \brief Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	inline grid_dist_key_dx<dim, typename device_grid::base_key> get() const
	{
		return grid_dist_key_dx<dim,typename device_grid::base_key>(g_c,a_it.get());
	}

	/*! \brief it return the stop point of the iterator
	 *
	 * The stop point of the iterator is just the grid size
	 *
	 * \return the stop point
	 *
	 */
	inline grid_key_dx<dim> getStop() const
	{
		return stop;
	}

	/*! \brief it return the start point of the iterator
	 *
	 * The start point of the iterator is the point with all coordinates zeros
	 *
	 * \return the start point
	 *
	 */
	inline grid_key_dx<dim> getStart() const
	{
		grid_key_dx<dim> start;

		start.zero();

		return start;
	}

	/*! \brief Get the boxes
	 *
	 *  Get the boxes that define the local grids
	 *
	 * \return Vector of local boxes
	 *
	 */
	inline const openfpm::vector<GBoxes<device_grid::dims>> & getGBoxes()
	{
		return gdb_ext;
	}

	/*! \brief Convert a g_dist_key_dx into a global key
	 *
	 * \see grid_dist_key_dx
	 * \see grid_dist_iterator
	 *
	 * \param k key position in local coordinates
	 *
	 * \return the global position in the grid
	 *
	 */
	inline grid_key_dx<dim> getGKey(const grid_dist_key_dx<dim,typename device_grid::base_key> & k)
	{
		// Get the sub-domain id
		size_t sub_id = k.getSub();

		grid_key_dx<dim> k_glob = k.getKey();

		// shift
		k_glob = k_glob + gdb_ext.get(sub_id).origin;

		return k_glob;
	}

	/*! \brief Return the stencil point offset
	 *
	 * \tparam id
	 *
	 * \return linearized distributed key
	 *
	 */
	template<unsigned int id> inline grid_dist_lin_dx getStencil()
	{
		return grid_dist_lin_dx(g_c,a_it.template getStencil<id>());
	}
};



#endif /* GRID_DIST_ID_ITERATOR_SUB_HPP_ */
