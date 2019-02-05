/*
 * grid_dist_id_iterator_util.hpp
 *
 *  Created on: Jan 6, 2017
 *      Author: i-bird
 */

#ifndef SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATOR_UTIL_HPP_
#define SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATOR_UTIL_HPP_


/*! \brief compute the subset where it has to iterate
 *
 * \param g_c Actual grid
 * \param start iterator start in global coordinate
 * \param stop iterator stop in global coordinate
 * \param start_c adjusted start point for the grid g_c
 * \param stop_c adjusted stop point for the grid g_c
 *
 * \return false if the sub-set does not contain points
 *
 */
template<typename Decomposition> static inline bool compute_subset(const openfpm::vector<GBoxes<Decomposition::dims>> & gdb_ext, size_t g_c, grid_key_dx<Decomposition::dims> & start, grid_key_dx<Decomposition::dims> & stop, grid_key_dx<Decomposition::dims> & start_c, grid_key_dx<Decomposition::dims> & stop_c)
{
	// Intersect the grid keys

	for (size_t i = 0 ; i < Decomposition::dims ; i++)
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



#endif /* SRC_GRID_ITERATORS_GRID_DIST_ID_ITERATOR_UTIL_HPP_ */
