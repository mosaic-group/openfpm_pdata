/*
 * grid_dist_util.hpp
 *
 *  Created on: Jan 28, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_GRID_DIST_UTIL_HPP_
#define SRC_GRID_GRID_DIST_UTIL_HPP_

#include "NN/CellList/CellDecomposer.hpp"
#include "Decomposition/common.hpp"

/*! \brief get cellDecomposer parameters
 *
 * \tparam dim dimensionality
 *
 * \param c_g get the parameters of the cell decomposer
 * \param g_sz global grid parameters
 *
 */
template<unsigned int dim> void getCellDecomposerPar(size_t (& c_g)[dim], const size_t (& g_sz)[dim], const size_t (& bc)[dim])
{
	for (size_t i = 0 ; i < dim ; i++)
	{
		if (bc[i] == NON_PERIODIC)
			c_g[i] = (g_sz[i]-1 > 0)?(g_sz[i]-1):1;
		else
			c_g[i] = g_sz[i];
	}
}


/*! \brief Create NON_PERIODIC data structure
 *
 * \tparam dim Dimensionality
 *
 * \return structure that define the non periodicity of the grid
 *
 */
template<unsigned int dim> periodicity<dim> create_non_periodic()
{
	periodicity<dim> p;

	for(size_t i = 0 ; i < dim ; i++)
		p.bc[i] = NON_PERIODIC;

	return p;
}

/*! \brief Create the gdb_ext
 *
 * It is a fundamental function, because it create the structure that store the information of the local grids. In
 * particular from the continuous decomposed domain it calculate the grid that each sub-domain store
 *
 * \param gdb_ext output Vector of Boxes that define the local grids extension
 * \param dec Decomposition
 * \param cd_sm CellDecomposer the size of cell is equal to the distance between grid points
 *
 */
template<int dim, typename Decomposition> inline void create_gdb_ext(openfpm::vector<GBoxes<Decomposition::dims>> & gdb_ext, Decomposition & dec, CellDecomposer_sm<Decomposition::dims,typename Decomposition::stype,shift<dim,typename Decomposition::stype>> & cd_sm)
{
	// Get the number of local grid needed
	size_t n_grid = dec.getNSubDomain();

	// Allocate the grids
	for (size_t i = 0 ; i < n_grid ; i++)
	{
		gdb_ext.add();

		// Get the local sub-domain (Grid conversion must be done with the domain P1 equivalent to 0.0)
		// consider that the sub-domain with point P1 equivalent to the domain P1 is a (0,0,0) in grid unit
		SpaceBox<Decomposition::dims, typename Decomposition::stype> sp = dec.getSubDomain(i);
		SpaceBox<Decomposition::dims, typename Decomposition::stype> sp_g = dec.getSubDomainWithGhost(i);

		// Because of round off we expand for safety the ghost area
		// std::nextafter return the next bigger or smaller representable floating
		// point number
		for (size_t i = 0 ; i < Decomposition::dims ; i++)
		{
			sp_g.setLow(i,std::nextafter(sp_g.getLow(i),sp_g.getLow(i) - 1.0));
			sp_g.setHigh(i,std::nextafter(sp_g.getHigh(i),sp_g.getHigh(i) + 1.0));
		}

		// Convert from SpaceBox<dim,St> to SpaceBox<dim,long int>
		SpaceBox<Decomposition::dims,long int> sp_t = cd_sm.convertDomainSpaceIntoGridUnits(sp,dec.periodicity());
		SpaceBox<Decomposition::dims,long int> sp_tg = cd_sm.convertDomainSpaceIntoGridUnits(sp_g,dec.periodicity());

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
 * \param sz Global grid grid size
 * \param domain Domain where the grid is defined
 * \param spacing Define the spacing of the grid
 * \param bc boundary conditions
 *
 */
template<int dim, typename Decomposition> inline void create_gdb_ext(openfpm::vector<GBoxes<dim>> & gdb_ext, Decomposition & dec, const size_t (& sz)[dim], const Box<Decomposition::dims,typename Decomposition::stype> & domain, typename Decomposition::stype (& spacing)[dim])
{
	// Create the cell decomposer
	CellDecomposer_sm<Decomposition::dims,typename Decomposition::stype, shift<Decomposition::dims,typename Decomposition::stype>> cd_sm;

	size_t cdp[dim];

	// Get the parameters to create a Cell-decomposer
	getCellDecomposerPar<Decomposition::dims>(cdp,sz,dec.periodicity());

	// Careful cd_sm require the number of cell
	cd_sm.setDimensions(domain,cdp,0);

	create_gdb_ext<dim,Decomposition>(gdb_ext,dec,cd_sm);

	// fill the spacing
	for (size_t i = 0 ; i < dim ; i++)
	{spacing[i] = cd_sm.getCellBox().getP2()[i];}
}

/*! \brief it store a box, its unique id and the sub-domain from where it come from
 *
 */
template<unsigned int dim> struct i_box_id
{
	//! Box
	::Box<dim,long int> box;

	//! id
	size_t g_id;

	//! r_sub id of the sub-domain in the sent list
	size_t r_sub;

	//! Sector where it live the linked external ghost box
	comb<dim> cmb;



	//! sub
	size_t sub;
};

/*! \brief it store an internal ghost box, the linked external ghost box and the sub-domain from where
 *  it come from as internal ghost box
 *
 */
template<unsigned int dim> struct i_lbox_id
{
	//! Box
	::Box<dim,long int> box;

	//! sub-domain id
	size_t sub;

	//! external ghost box linked to this internal ghost box
	size_t k;

	//! combination
	comb<dim> cmb;
};

/*! \brief It store the information about the external ghost box
 *
 *
 */
template <unsigned int dim> struct e_box_id
{
	//! Box defining the external ghost box in global coordinates
	::Box<dim,long int> g_e_box;

	//! Box defining the external ghost box in local coordinates
	::Box<dim,long int> l_e_box;

	//! Sector position of the external ghost
	comb<dim> cmb;

	//! Id
	size_t g_id;

	//! sub_id in which sub-domain this box live
	size_t sub;
};

/*! \brief It store the information about the local external ghost box
 *
 *
 */
template <unsigned int dim> struct e_lbox_id
{
	//! Box defining the external ghost box in local coordinates
	::Box<dim,long int> box;

	//! Has this external ghost box initialized
	bool initialized = false;

	//! sub_id in which sub-domain this box live
	size_t sub;

	//! external ghost box linked to this internal ghost box
	size_t k;

	//! Sector position of the local external ghost box
	comb<dim> cmb;
};

/*! \brief Per-processor Internal ghost box
 *
 */
template <unsigned int dim> struct ip_box_grid
{
	// ghost in grid units
	openfpm::vector<i_box_id<dim>> bid;

	//! processor id
	size_t prc;
};

/*! \brief local Internal ghost box
 *
 */
template <unsigned int dim> struct i_lbox_grid
{
	// ghost in grid units
	openfpm::vector<i_lbox_id<dim>> bid;
};

/*! \brief Per-processor external ghost box
 *
 */
template <unsigned int dim>struct ep_box_grid
{
	// ghost in grid units
	openfpm::vector<e_box_id<dim>> bid;

	//! processor id
	size_t prc;
};

/*! \brief Per-processor external ghost box
 *
 */
template <unsigned int dim> struct e_lbox_grid
{
	// ghost in grid units
	openfpm::vector<e_lbox_id<dim>> bid;
};

#endif /* SRC_GRID_GRID_DIST_UTIL_HPP_ */
