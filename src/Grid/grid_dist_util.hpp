/*
 * grid_dist_util.hpp
 *
 *  Created on: Jan 28, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_UTIL_HPP_
#define SRC_GRID_GRID_DIST_UTIL_HPP_

#include "NN/CellList/CellDecomposer.hpp"

/*! \brief Create the gdb_ext
 *
 * \param gdb_ext Vector of Boxes that define the local grids extension
 * \param dec Decomposition
 * \param cd_sm CellDecomposer the size of cell is equal to the distance between grid points
 *
 */
template<int dim, typename Decomposition> inline void create_gdb_ext(openfpm::vector<GBoxes<Decomposition::dims>> & gdb_ext, Decomposition & dec, CellDecomposer_sm<Decomposition::dims,typename Decomposition::stype> & cd_sm)
{
	Box<Decomposition::dims, typename Decomposition::stype> g_rnd_box;
	for (size_t i = 0 ; i < Decomposition::dims ; i++)	{g_rnd_box.setHigh(i,0.5); g_rnd_box.setLow(i,-0.5);}

	// Get the number of local grid needed
	size_t n_grid = dec.getNLocalHyperCube();

	// Allocate the grids
	for (size_t i = 0 ; i < n_grid ; i++)
	{
		gdb_ext.add();

		// Get the local hyper-cube
		SpaceBox<Decomposition::dims, typename Decomposition::stype> sp = dec.getLocalHyperCube(i);
		SpaceBox<Decomposition::dims, typename Decomposition::stype> sp_g = dec.getSubDomainWithGhost(i);

		// Convert from SpaceBox<dim,St> to SpaceBox<dim,long int>
		SpaceBox<Decomposition::dims,long int> sp_t = cd_sm.convertDomainSpaceIntoGridUnits(sp);
		SpaceBox<Decomposition::dims,long int> sp_tg = cd_sm.convertDomainSpaceIntoGridUnits(sp_g);

		//! Save the origin of the sub-domain of the local grid
		gdb_ext.last().origin = sp_tg.getP1();

		// save information about the local grid: domain box seen inside the domain + ghost box (see GDBoxes for a visual meaning)
		// and where the GDBox start, or the origin of the local grid (+ghost) in global coordinate
		gdb_ext.last().Dbox = sp_t;
		gdb_ext.last().Dbox -= sp_tg.getP1();

		gdb_ext.last().GDbox = sp_tg;
		gdb_ext.last().GDbox -= sp_tg.getP1();
	}
}

/*! \brief Create the gdb_ext
 *
 * \param gdb_ext Vector of Boxes that define the local grids extension
 * \param dec Decomposition
 * \param Global grid grid size
 *
 */
template<int dim, typename Decomposition> inline void create_gdb_ext(openfpm::vector<GBoxes<dim>> & gdb_ext, Decomposition & dec, const size_t (& sz)[dim], const Box<Decomposition::dims,typename Decomposition::stype> & domain)
{
	// Create the cell decomposer

	CellDecomposer_sm<Decomposition::dims,typename Decomposition::stype> cd_sm;
	cd_sm.setDimensions(domain,sz,0);

	create_gdb_ext<dim,Decomposition>(gdb_ext,dec,cd_sm);
}

#endif /* SRC_GRID_GRID_DIST_UTIL_HPP_ */
