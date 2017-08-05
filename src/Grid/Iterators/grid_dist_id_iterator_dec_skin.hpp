/*
 * grid_dist_id_iterator_dec_skin.hpp
 *
 *  Created on: Jan 4, 2017
 *      Author: i-bird
 */

#ifndef SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATOR_DEC_SKIN_HPP_
#define SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATOR_DEC_SKIN_HPP_


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
template<typename Decomposition>
class grid_dist_id_iterator_dec_skin : protected grid_skin_iterator_bc<Decomposition::dims>
{
	//! a_its element in this moment selected
	size_t a_its_p;

	//! grid list counter
	size_t g_c;

	//! Extension of each grid: domain and ghost + domain
	openfpm::vector<GBoxes<Decomposition::dims>> gdb_ext;

	//! Internal grid sub-iterator
	grid_key_dx_iterator_sub<Decomposition::dims> a_it;

	//! Internal struct
	struct gp_sub
	{
		//! from which grid this iterator come from
		size_t gc;

		//! Iterator
		grid_key_dx_iterator_sub<Decomposition::dims> it;

		/*! \brief constructor
		 *
		 * \param gc sub-domain
		 * \param it iterator
		 *
		 */
		gp_sub(size_t gc, grid_key_dx_iterator_sub<Decomposition::dims> && it)
		:gc(gc),it(it)
		{}
	};

	//! Actual sub-iterators
	openfpm::vector<gp_sub> a_its;

	//! Spacing
	typename Decomposition::stype spacing[Decomposition::dims];


	/*! \brief from g_c increment g_c until you find a valid grid
	 *
	 */
	void selectValidGrid()
	{
		if (a_its_p < a_its.size())
		{
			g_c = a_its.get(a_its_p).gc;

			a_it.reinitialize(a_its.get(a_its_p).it);
		}
		else
			g_c = gdb_ext.size();
	}

	/*! \brief construct sub-iterators
	 *
	 *
	 *
	 */
	void construct_sub_it()
	{
		// Construct the sub iterators
		for (size_t i = 0 ; i < 2*Decomposition::dims; i++)
		{
			for (size_t gc = 0 ; gc < gdb_ext.size() ; gc++)
			{
				grid_key_dx<Decomposition::dims> start = this->sub_it[i].getStart();
				grid_key_dx<Decomposition::dims> stop = this->sub_it[i].getStop();

				grid_key_dx<Decomposition::dims> start_c;
				grid_key_dx<Decomposition::dims> stop_c;

				if (compute_subset<Decomposition>(gdb_ext,gc,start,stop,start_c,stop_c) == true)
				{
					// Convert global coordinate start_c stop_c into local
					// and calculate the grid sizes
					size_t sz[Decomposition::dims];
					for (size_t j = 0 ; j < Decomposition::dims ; j++)
						sz[j] = gdb_ext.get(gc).GDbox.getHigh(j) + 1;

					grid_sm<Decomposition::dims,void> g_sm(sz);

					// Non empty sub-set
					a_its.add(gp_sub(gc,grid_key_dx_iterator_sub<Decomposition::dims>(g_sm,start_c,stop_c)));
				}
			}
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


	/*! \brief Copy constructor
	*
	* \param tmp iterator to copy
	*
	*/
	grid_dist_id_iterator_dec_skin(const grid_dist_id_iterator_dec_skin<Decomposition> & tmp)
	:grid_skin_iterator_bc<Decomposition::dims>(tmp),a_its_p(0)
	{
		this->operator=(tmp);
	}

	/*! \brief Copy constructor
	*
	* \param tmp iterator to copy
	*
	*/
	grid_dist_id_iterator_dec_skin(grid_dist_id_iterator_dec_skin<Decomposition> && tmp)
	:grid_skin_iterator_bc<Decomposition::dims>(tmp),a_its_p(0)
	{
		this->operator=(tmp);
	}

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param dec Decomposition
	 * \param g_sm grid size on each direction
	 * \param A box A (must be contained into B)
	 * \param B box B
	 * \param bc boundary conditions
	 *
	 */
	grid_dist_id_iterator_dec_skin(Decomposition & dec,
			                       const grid_sm<Decomposition::dims,void> & g_sm,
								   const Box<Decomposition::dims,size_t> & A,
								   const Box<Decomposition::dims,size_t> & B,
								   const size_t (& bc)[Decomposition::dims])
	:grid_skin_iterator_bc<Decomposition::dims>(g_sm,A,B,bc),a_its_p(0),g_c(0)
	{
		// From the decomposition construct gdb_ext
		create_gdb_ext<Decomposition::dims,Decomposition>(gdb_ext,dec,g_sm.getSize(),dec.getDomain(),spacing);

		// This iterato only work if A is contained into B
		if (A.isContained(B) == false)
			std::cout << __FILE__ << ":" << __LINE__ << ", Error Box A must be contained into box B" << std::endl;

		// construct sub_iterators
		construct_sub_it();

		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}

	//! Destructor
	~grid_dist_id_iterator_dec_skin()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return itself
	 *
	 */
	inline grid_dist_id_iterator_dec_skin<Decomposition> & operator++()
	{
		++a_it;

		// check if a_it is at the end

		if (a_it.isNext() == true)
			return *this;
		else
		{
			// switch to the new grid
			a_its_p++;

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
	 * \param i dimension
	 *
	 * \return the spacing
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

		if (k_glob.get(0) > 11)
		{
			int debug = 0;
			debug++;
		}

		return k_glob;
	}

	/*! \brief Copy operator=
	*
	* \param tmp iterator to copy
	*
	* \return itself
	*
	*/
	grid_dist_id_iterator_dec_skin<Decomposition> & operator=(const grid_dist_id_iterator_dec_skin<Decomposition> & tmp)
	{
		a_its_p = tmp.a_its_p;
		g_c = tmp.g_c;
		gdb_ext = tmp.gdb_ext;
		a_its = tmp.a_its;

		for (size_t i = 0 ; i < Decomposition::dims ; i++)
			spacing[i] = tmp.spacing[i];

		a_it.reinitialize(tmp.a_it);

		return *this;
	}

	/*! \brief Copy operator=
	*
	* \param tmp iterator to copy
	*
	* \return itself
	*
	*/
	grid_dist_id_iterator_dec_skin<Decomposition> & operator=(grid_dist_id_iterator_dec_skin<Decomposition> && tmp)
	{
		a_its_p = tmp.a_its_p;
		g_c = tmp.g_c;
		gdb_ext = tmp.gdb_ext;
		a_its = tmp.a_its;

		for (size_t i = 0 ; i < Decomposition::dims ; i++)
			spacing[i] = tmp.spacing[i];

		a_it.reinitialize(tmp.a_it);

		return *this;
	}
};


#endif /* SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATOR_DEC_SKIN_HPP_ */
