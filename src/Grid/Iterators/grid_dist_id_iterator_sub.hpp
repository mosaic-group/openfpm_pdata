/*
 * grid_dist_id_iterator_sub.hpp
 *
 *  Created on: Oct 14, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_ID_ITERATOR_SUB_HPP_
#define SRC_GRID_GRID_DIST_ID_ITERATOR_SUB_HPP_


/*! \brief Distributed grid iterator
 *
 * Iterator across the local elements of the distributed grid
 *
 * \tparam dim dimensionality of the grid
 * \tparam device_grid type of basic grid
 * \tparam impl implementation
 *
 */
template<unsigned int dim, typename device_grid>
class grid_dist_iterator_sub
{
	// sub_set of the grid where to iterate
	struct sub_set
	{
		//! start point where iterate
		grid_key_dx<dim> start;
		//! stop point where iterate
		grid_key_dx<dim> stop;
	};

	//! grid list counter
	size_t g_c;

	//! List of the grids we are going to iterate
	const openfpm::vector<device_grid> & gList;

	//! Extension of each grid: domain and ghost + domain
	const openfpm::vector_ofp<GBoxes<device_grid::dims>> & gdb_ext;

	//! Actual iterator
	grid_key_dx_iterator_sub<dim> a_it;

	//! start key
	grid_key_dx<dim> start;

	//! stop key
	grid_key_dx<dim> stop;

	/*! \brief compute the subset where it has to iterate
	 *
	 * \param gc Actual grid
	 * \param start_c adjusted start point for the grid g_c
	 * \param stop_c adjusted stop point for the grid g_c
	 *
	 * \return false if the sub-set does not contain points
	 *
	 */
	bool compute_subset(size_t gc, grid_key_dx<dim> & start_c, grid_key_dx<dim> & stop_c)
	{
		// Intersect the grid keys

		for (size_t i = 0 ; i < dim ; i++)
		{
			long int start_p = gdb_ext.get(g_c).Dbox.getP1().get(i) + gdb_ext.get(g_c).origin.get(i);
			long int stop_p = gdb_ext.get(g_c).Dbox.getP2().get(i) + gdb_ext.get(g_c).origin.get(i);
			if (start.get(i) <= start_p)
				start_c.set_d(i,gdb_ext.get(g_c).Dbox.getP1().get(i));
			else if (start.get(i) <= stop_p)
				start_c.set_d(i,start.get(i) - gdb_ext.get(g_c).origin.get(i));
			else
				return false;

			if (stop.get(i) >= stop_p)
				stop_c.set_d(i,gdb_ext.get(g_c).Dbox.getP2().get(i));
			else if (stop.get(i) >= start_p)
				stop_c.set_d(i,stop.get(i) - gdb_ext.get(g_c).origin.get(i));
			else
				return false;
		}

		return true;
	}

	/*! \brief from g_c increment g_c until you find a valid grid
	 *
	 */
	void selectValidGrid()
	{
		// start and stop for the subset grid
		grid_key_dx<dim> start_c;
		grid_key_dx<dim> stop_c;

		// When the grid has size 0 potentially all the other informations are garbage
		while (g_c < gList.size() &&
			   (gList.get(g_c).size() == 0 || gdb_ext.get(g_c).Dbox.isValid() == false || compute_subset(g_c,start_c,stop_c) == false ))
		{g_c++;}

		// get the next grid iterator
		if (g_c < gList.size())
		{
			a_it.reinitialize(gList.get(g_c).getIterator(start_c,stop_c));
		}
	}

	public:

	/*! \brief Copy operator=
	*
	* \param tmp iterator to copy
	*
	*/
	grid_dist_iterator_sub<dim,device_grid> & operator=(const grid_dist_iterator_sub<dim,device_grid> & tmp)
	{
		g_c = tmp.g_c;
		gList = tmp.gList;
		gdb_ext = tmp.gdb_ext;
		start = tmp.start;
		stop = tmp.stop;
		a_it.reinitialize(tmp.a_it);

		return *this;
	}

	/*! \brief Copy constructor
	*
	* \param tmp iterator to copy
	*
	*/
	grid_dist_iterator_sub(const grid_dist_iterator_sub<dim,device_grid> & tmp)
	:g_c(tmp.g_c),gList(tmp.gList),gdb_ext(tmp.gdb_ext),start(tmp.start),stop(tmp.stop)
	{
		// get the next grid iterator
		if (g_c < gList.size())
		{
			a_it.reinitialize(tmp.a_it);
		}
	}

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param start position
	 * \param stop position
	 * \param gk std::vector of the local grid
	 * \param gdb_ext information about the local grids
	 *
	 */
	grid_dist_iterator_sub(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop ,const openfpm::vector<device_grid> & gk, const openfpm::vector_ofp<GBoxes<device_grid::dims>> & gdb_ext)
	:g_c(0),gList(gk),gdb_ext(gdb_ext),start(start),stop(stop)
	{
		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}

	// Destructor
	~grid_dist_iterator_sub()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */

	inline grid_dist_iterator_sub<dim,device_grid> & operator++()
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
	inline bool isNext()
	{
		// If there are no other grid stop

		if (g_c >= gList.size())
			return false;

		return true;
	}

	/*! \brief Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	inline grid_dist_key_dx<dim> get()
	{
		return grid_dist_key_dx<dim>(g_c,a_it.get());
	}

	/*! \brief Convert a g_dist_key_dx into a global key
	 *
	 * \see grid_dist_key_dx
	 * \see grid_dist_iterator
	 *
	 * \return the global position in the grid
	 *
	 */
	inline grid_key_dx<dim> getGKey(const grid_dist_key_dx<dim> & k)
	{
		// Get the sub-domain id
		size_t sub_id = k.getSub();

		grid_key_dx<dim> k_glob = k.getKey();

		// shift
		k_glob = k_glob + gdb_ext.get(sub_id).origin;

		return k_glob;
	}

	/* \brief Get the starting point of the grid iterator
	 *
	 * \return the starting point
	 *
	 */
	inline grid_key_dx<dim> getStart() const
	{
		return start;
	}

	/* \brief Get the stop point of the grid iterator
	 *
	 * \return the stop point
	 *
	 */
	inline grid_key_dx<dim> getStop() const
	{
		return stop;
	}

	/*! \brief Return the number of local grids
	 *
	 *
	 */
	inline size_t N_loc_grid()
	{
		return gList.size();
	}

	/*! \brief Return the component j of the starting point (P1) of the domain part
	 *         for the local grid i
	 *
	 * \param i local grid
	 * \param j dimension
	 *
	 *
	 */
	inline size_t loc_grid_info_start(size_t i,size_t j)
	{
		return gdb_ext.get(i).DBox.getLow(i);
	}

	/*! \brief Return the component j of the stop point (P2) of the domain part
	 *         for the local grid i
	 *
	 * \param i local grid
	 * \param j dimension
	 *
	 *
	 */
	inline size_t loc_grid_info_size(size_t i,size_t j)
	{
		return gdb_ext.get(i).GDBox.getHigh(i);
	}
};


//////// MACRO in 3D

#define WHILE_M(grid,stencil) auto & ginfo = grid.getLocalGridsInfo();\
								 for (size_t s = 0 ; s < grid.getN_loc_grid() ; s++)\
								 {\
									 auto it = grid.get_loc_grid_iterator_stencil(s,stencil);\
\
										int lo[3] = {(int)ginfo.get(s).Dbox.getLow(0),(int)ginfo.get(s).Dbox.getLow(1),(int)ginfo.get(s).Dbox.getLow(2)};\
										int hi[3] = {(int)ginfo.get(s).Dbox.getHigh(0),(int)ginfo.get(s).Dbox.getHigh(1),(int)ginfo.get(s).Dbox.getHigh(2)};\
\
										int uhi[3] = {(int)ginfo.get(s).GDbox.getHigh(0),(int)ginfo.get(s).GDbox.getHigh(1),(int)ginfo.get(s).GDbox.getHigh(2)};\
\
										int sx = uhi[0]+1;\
										int sxsy = (uhi[0]+1)*(uhi[1]+1);

#define ITERATE_3D_M(n_pt)			int i = lo[2];\
									for ( ; i <= hi[2] ; i+=1)\
									{\
										int j = lo[1];\
										for ( ; j <= hi[1] ; j+=1)\
										{\
											int k = lo[0];\
											for ( ; k <= hi[0] ; k+=n_pt)\
											{


#define GET_GRID_M(grid)	grid.get_loc_grid(s);


#define END_LOOP_M(n_pt) 					it.private_sum<n_pt>();\
								}\
								it.private_adjust( - k + sx + lo[0]);\
							}\
							it.private_adjust(- j*sx + sxsy + lo[1]*sx);\
						}\
					}

#endif /* SRC_GRID_GRID_DIST_ID_ITERATOR_SUB_HPP_ */
