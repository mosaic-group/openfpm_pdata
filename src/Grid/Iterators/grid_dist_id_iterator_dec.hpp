/*
 * grid_dist_id_iterator_dec.hpp
 *
 *  Created on: Jan 27, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_ID_ITERATOR_DEC_HPP_
#define SRC_GRID_GRID_DIST_ID_ITERATOR_DEC_HPP_

#include "grid_dist_id_iterator.hpp"
#include "Grid/grid_dist_util.hpp"
#include "grid_dist_id_iterator_util.hpp"

/*! \brief Given the decomposition it create an iterator
 *
 * Iterator across the local elements of the distributed grid
 *
 * \tparam dec Decomposition type
 *
 */
template<typename Decomposition, bool ghost_or_domain = false>
class grid_dist_id_iterator_dec
{
	//! grid list counter
	size_t g_c;

	//! Extension of each grid: domain and ghost + domain
	openfpm::vector<GBoxes<Decomposition::dims>> gdb_ext;

	//! Actual iterator
	grid_key_dx_iterator_sub<Decomposition::dims> a_it;

	//! start key
	grid_key_dx<Decomposition::dims> start;

	//! stop key
	grid_key_dx<Decomposition::dims> stop;

	//! Spacing
	typename Decomposition::stype spacing[Decomposition::dims];

	//! Domain
	Box<Decomposition::dims,typename Decomposition::stype> domain;

	/*! \brief from g_c increment g_c until you find a valid grid
	 *
	 */
	void selectValidGrid()
	{
		// start and stop for the subset grid
		grid_key_dx<Decomposition::dims> start_c;
		grid_key_dx<Decomposition::dims> stop_c;

		// When the grid has size 0 potentially all the other informations are garbage
		while (g_c < gdb_ext.size() &&
			   (gdb_ext.get(g_c).Dbox.isValid() == false || compute_subset<Decomposition,ghost_or_domain>(gdb_ext,g_c,start,stop,start_c,stop_c) == false ))
		{g_c++;}

		// get the next grid iterator
		if (g_c < gdb_ext.size())
		{
			// Calculate the resolution of the local grid
			size_t sz[Decomposition::dims];
			for (size_t i = 0 ; i < Decomposition::dims ; i++)
				sz[i] = gdb_ext.get(g_c).GDbox.getP2()[i] + 1;

			grid_sm<Decomposition::dims,void> g_sm(sz);
			a_it.reinitialize(grid_key_dx_iterator_sub<Decomposition::dims>(g_sm,start_c,stop_c));
		}
	}

	/*! \brief Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	inline grid_dist_key_dx<Decomposition::dims> get_int()
	{
		return grid_dist_key_dx<Decomposition::dims>(g_c,a_it.get());
	}

	public:

	/*! \brief Copy operator=
	*
	* \param tmp iterator to copy
	*
	* \return itself
	*
	*/
	grid_dist_id_iterator_dec<Decomposition> & operator=(const grid_dist_id_iterator_dec<Decomposition> & tmp)
	{
		g_c = tmp.g_c;
		gdb_ext = tmp.gdb_ext;
		a_it.reinitialize(tmp.a_it);

		start = tmp.start;
		stop = tmp.stop;

		domain = tmp.domain;

		return *this;
	}

	/*! \brief Copy constructor
	*
	* \param tmp iterator to copy
	*
	*/
	grid_dist_id_iterator_dec(const grid_dist_id_iterator_dec<Decomposition> & tmp)
	{
		this->operator=(tmp);
	}

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param dec Decomposition
	 * \param sz size of the grid
	 *
	 */
	grid_dist_id_iterator_dec(Decomposition & dec, const size_t (& sz)[Decomposition::dims])
	:g_c(0)
	{
		domain = dec.getDomain();

		// Initialize start and stop
		start.zero();
		for (size_t i = 0 ; i < Decomposition::dims ; i++)
			stop.set_d(i,sz[i]-1);

		// From the decomposition construct gdb_ext
		create_gdb_ext<Decomposition::dims,Decomposition>(gdb_ext,dec,sz,dec.getDomain(),spacing);

		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param dec Decomposition
	 * \param sz size of the grid
	 * \param start point
	 * \param stop point
	 *
	 */
	grid_dist_id_iterator_dec(Decomposition & dec, const size_t (& sz)[Decomposition::dims], grid_key_dx<Decomposition::dims> start, grid_key_dx<Decomposition::dims> stop)
	:g_c(0),start(start),stop(stop)
	{
		domain = dec.getDomain();

		// From the decomposition construct gdb_ext
		create_gdb_ext<Decomposition::dims,Decomposition>(gdb_ext,dec,sz,dec.getDomain(),spacing);

		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}

	// Destructor
	~grid_dist_id_iterator_dec()
	{
	}

	/*! \brief Return true if we point to a valid grid
	 *
	 * \return true if valid grid
	 *
	 */
	inline bool isNextGrid()
	{
		return g_c < gdb_ext.size();
	}

	/*! \brief Return the index of the grid in which we are iterating
	 *
	 *
	 */
	inline size_t getGridId()
	{
		return g_c;
	}

	/*! \brief next grid
	 *
	 *
	 */
	inline void nextGrid()
	{
		g_c++;
		selectValidGrid();
	}

	/*! \brief Return the actual pointed grid
	 *
	 * \return the grid index
	 *
	 */
	inline Box<Decomposition::dims,size_t> getGridBox()
	{
		Box<Decomposition::dims,size_t> bx;

		auto start = a_it.getStart();
		auto stop = a_it.getStop();

		for (int i = 0 ; i < Decomposition::dims ; i++)
		{
			bx.setHigh(i,stop.get(i));
			bx.setLow(i,start.get(i));
		}

		return bx;
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	inline grid_dist_id_iterator_dec<Decomposition,ghost_or_domain> & operator++()
	{
		++a_it;

		// check if a_it is at the end

		if (a_it.isNext() == true)
		{return *this;}
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
	inline bool isNext()
	{
		// If there are no other grid stop

		if (g_c >= gdb_ext.size())
			return false;

		return true;
	}

	/*! \brief Get the spacing of the grid
	 *
	 * \param i
	 *
	 */
	inline typename Decomposition::stype getSpacing(size_t i)
	{
		return spacing[i];
	}

	/*! \brief Get the actual global key of the grid
	 *
	 *
	 * \return the global position in the grid
	 *
	 */
	inline grid_key_dx<Decomposition::dims> get()
	{
		const grid_dist_key_dx<Decomposition::dims> k = get_int();

		// Get the sub-domain id
		size_t sub_id = k.getSub();

		grid_key_dx<Decomposition::dims> k_glob = k.getKey();

		// shift
		k_glob = k_glob + gdb_ext.get(sub_id).origin;

		return k_glob;
	}

	/*! \brief Return the point coordinates
	 *
	 * \return the point
	 *
	 */
	inline Point<Decomposition::dims,typename Decomposition::stype> getPoint()
	{
		Point<Decomposition::dims,typename Decomposition::stype> p;
		auto key = this->get();

		for (int i = 0 ; i < Decomposition::dims ; i++)
		{
			p.get(i) = spacing[i] * key.get(i) + domain.getLow(i);
		}

		return p;
	}

	/*! \brief Get the actual grid key for a distributed grid
	 *
	 * \note if you are using this iterator and you need the position for grid_dist_id use this
	 *       function to retrieve the grid point
	 *
	 * \return the grid point
	 *
	 */
	inline grid_dist_key_dx<Decomposition::dims> get_dist()
	{
		return get_int();
	}

	/*! \brief Get the starting point of the sub-grid we are iterating
	 *
	 * \return the starting point
	 *
	 */
	inline grid_key_dx<Decomposition::dims> getStart()
	{
		return start;
	}

	/*! \brief Get the starting point of the sub-grid we are iterating
	 *
	 * \return the stop point
	 *
	 */
	inline grid_key_dx<Decomposition::dims> getStop()
	{
		return stop;
	}
};


#endif /* SRC_GRID_GRID_DIST_ID_ITERATOR_DEC_HPP_ */
