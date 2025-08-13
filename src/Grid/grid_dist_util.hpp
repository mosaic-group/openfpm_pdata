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

template<unsigned int dim>
size_t get_gdb_ext(const openfpm::vector<GBoxes<dim>> & gdb_ext, size_t start , size_t stop, size_t k)
{
	for (size_t i = start ; i < stop ; i++)
	{
		if (gdb_ext.get(i).k == k)
		{return i;}
	}

	// not found
	return (size_t)-1;
}

/*! \brief Add a box to gdb_ext
 *
 * \param on which gdb_ext array add the box
 * \param k from which definition box it come from
 *        This parameter make sense only when the grid is not defined everywhere
 * \param grid of the sub-domain
 * \param grid of the sub-domain inclusive of ghost
 *
 */
template<unsigned int dim>
void add_to_gdb_ext(openfpm::vector<GBoxes<dim>> & gdb_ext, size_t k, Box<dim,long int> & sp_t, Box<dim,long int> & sp_tg)
{
	// Add gdb_ext
	gdb_ext.add();

	//! Save the origin of the sub-domain of the local grid
	gdb_ext.last().origin = sp_tg.getP1();

	// save information about the local grid: domain box seen inside the domain + ghost box (see GDBoxes for a visual meaning)
	// and where the GDBox start, or the origin of the local grid (+ghost) in global coordinate
	gdb_ext.last().Dbox = sp_t;
	gdb_ext.last().Dbox -= sp_tg.getP1();

	gdb_ext.last().GDbox = sp_tg;
	gdb_ext.last().GDbox -= sp_tg.getP1();

	gdb_ext.last().k = k;
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
template<int dim, typename Decomposition>
inline void create_gdb_ext(openfpm::vector<GBoxes<Decomposition::dims>> & gdb_ext,
						   openfpm::vector<size_t> & gdb_ext_markers,
		                   Decomposition & dec,
						   CellDecomposer_sm<Decomposition::dims,typename Decomposition::stype,shift<dim,typename Decomposition::stype>> & cd_sm,
						   openfpm::vector<Box<Decomposition::dims,long int>> & bx_create,
						   const Ghost<Decomposition::dims,long int> & exp,
						   bool use_bx_def)
{
	gdb_ext.clear();
	gdb_ext_markers.clear();

	// Get the number of local grid needed
	size_t n_grid = dec.getNSubDomain();

	// Allocate the grids
	for (size_t i = 0 ; i < n_grid ; i++)
	{
		gdb_ext_markers.add(gdb_ext.size());

		// Get the local sub-domain (Grid conversion must be done with the domain P1 equivalent to 0.0)
		// consider that the sub-domain with point P1 equivalent to the domain P1 is a (0,0,0) in grid unit
		Box<Decomposition::dims, typename Decomposition::stype> sp = dec.getSubDomain(i);
		Box<Decomposition::dims, typename Decomposition::stype> sp_g = dec.getSubDomainWithGhost(i);

		// Because of round off we expand for safety the ghost area
		// std::nextafter return the next bigger or smaller representable floating
		// point number
		for (size_t i = 0 ; i < Decomposition::dims ; i++)
		{
			sp_g.setLow(i,std::nextafter(sp_g.getLow(i),sp_g.getLow(i) - 1.0));
			sp_g.setHigh(i,std::nextafter(sp_g.getHigh(i),sp_g.getHigh(i) + 1.0));
		}

		// Convert from Box<dim,St> to Box<dim,long int>
		Box<Decomposition::dims,long int> sp_t = cd_sm.convertDomainSpaceIntoGridUnits(sp,dec.periodicity());
		Box<Decomposition::dims,long int> sp_tg = cd_sm.convertDomainSpaceIntoGridUnits(sp_g,dec.periodicity());

		for (size_t i = 0 ; i < Decomposition::dims ; i++)
		{
			if (sp_t.getLow(i) < sp_tg.getLow(i))
			{sp_tg.setLow(i,sp_t.getLow(i));}
			if (sp_t.getHigh(i) > sp_tg.getHigh(i))
			{sp_tg.setHigh(i,sp_t.getHigh(i));}
		}

		if (use_bx_def == true)
		{
			// intersect the sub-domain with all the boxes

			for (size_t k = 0 ; k < bx_create.size() ; k++)
			{
				Box<Decomposition::dims, long int> inte;

				if (sp_t.Intersect(bx_create.get(k),inte) == true)
				{
					// Ok we have a sub-domain now we have to create the ghost part.
					// The ghost part is created converting the bx_def into a continuous
					// box expanding this box by the ghost and intersecting this box
					// with the sub-domain. This is one way to get a ghost area consistent
					// with the construction of the external and internal ghost boxes,
					// always calculated in continuous from the decomposition.
					//

					Box<Decomposition::dims,typename Decomposition::stype> output;
					Box<Decomposition::dims,typename Decomposition::stype> bx_wg = cd_sm.convertCellUnitsIntoDomainSpaceMiddle(bx_create.get(k));
					bx_wg.enlarge(dec.getGhost());
					bx_wg.Intersect(sp_g,output);

					Box<Decomposition::dims,long int> sp_t2 = inte;
					Box<Decomposition::dims,long int> sp_tg2 =  cd_sm.convertDomainSpaceIntoGridUnits(output,dec.periodicity());

					add_to_gdb_ext(gdb_ext,k,sp_t2,sp_tg2);
				}
			}
		}
		else
		{
			add_to_gdb_ext(gdb_ext,0,sp_t,sp_tg);
		}
	}

	gdb_ext_markers.add(gdb_ext.size());
}

/*! \brief Create the gdb_ext
 *
 * \param gdb_ext Vector of Boxes that define the local grids extension
 * \param gdb_ext_markers filled with sub-domain markers
 *                see gdb_ext_markers in grid_dist_id for an explanation
 * \param dec Decomposition
 * \param sz Global grid grid size
 * \param domain Domain where the grid is defined
 * \param spacing Define the spacing of the grid
 * \param bc boundary conditions
 *
 */
template<int dim, typename Decomposition>
inline void create_gdb_ext(openfpm::vector<GBoxes<dim>> & gdb_ext,
		                   Decomposition & dec,
						   const size_t (& sz)[dim],
						   const Box<Decomposition::dims,typename Decomposition::stype> & domain,
						   typename Decomposition::stype (& spacing)[dim])
{
	// Create the cell decomposer
	CellDecomposer_sm<Decomposition::dims,typename Decomposition::stype, shift<Decomposition::dims,typename Decomposition::stype>> cd_sm;

	size_t cdp[dim];

	// Get the parameters to create a Cell-decomposer
	getCellDecomposerPar<Decomposition::dims>(cdp,sz,dec.periodicity());

	// Careful cd_sm require the number of cell
	cd_sm.setDimensions(domain,cdp,0);

	// create an empty vector of boxes
	openfpm::vector<Box<Decomposition::dims,long int>> empty;
	Ghost<Decomposition::dims,long int> zero(0);

	//! We are not interested on the markers
	openfpm::vector<size_t> unused;
	create_gdb_ext<dim,Decomposition>(gdb_ext,unused,dec,cd_sm,empty,zero,false);

	// fill the spacing
	for (size_t i = 0 ; i < dim ; i++)
	{spacing[i] = cd_sm.getCellBox().getP2()[i];}
}

/*! \brief It store the information about the external ghost box
 *
 *
 */
template <unsigned int dim> struct e_box_id
{
	//! Box defining the external ghost box in global coordinates
	::Box<dim,long int> g_e_box;

	//! Box defining the external ghost box in local coordinates for gdb_ext
	::Box<dim,long int> l_e_box;

	//! Box defining the external box in local coordinates for received box
	::Box<dim,long int> lr_e_box;

	//! Sector position of the external ghost
	comb<dim> cmb;

	//! Id
	size_t g_id;

	//! sub_id in which sub-domain this box live
	size_t sub;
};

/*! \brief flip box just convert and internal ghost box into an external ghost box and the other way around
 *
 * \param box to convert
 * \param cmb sector position of the box
 *
 * \return the converted box
 *
 */
template<unsigned int dim, typename T, typename idT>
Box<dim,long int> flip_box(const Box<dim,idT> & box, const comb<dim> & cmb, const grid_sm<dim,T> & ginfo)
{
	Box<dim,long int> flp;

	for (size_t i = 0 ; i < dim; i++)
	{
		if (cmb[i] == 0)
		{
			flp.setLow(i,box.getLow(i));
			flp.setHigh(i,box.getHigh(i));
		}
		else if (cmb[i] == 1)
		{
			flp.setLow(i,box.getLow(i) + ginfo.size(i));
			flp.setHigh(i,box.getHigh(i) + ginfo.size(i));
		}
		else if (cmb[i] == -1)
		{
			flp.setLow(i,box.getLow(i) - ginfo.size(i));
			flp.setHigh(i,box.getHigh(i) - ginfo.size(i));
		}
	}

	return flp;
}

/*! \brief convert to sub-domain id
 *
 * In case the grid is not defined everywhere the ids returned by getProcessorIGhostSub
 * or any function that return a sub-domain id does not match the ids in gdb_ext. This
 * function convert it to the correct one
 *
 * \param k sub-domain id to convert
 * \param def_id id of the box that define the real allocated grid
 * \param gdb_ext_markers markers for sub-domain id gdb_ext
 *
 */
template<unsigned int dim>
inline size_t convert_to_gdb_ext(size_t sub_id,
		                         size_t def_id,
								 openfpm::vector<GBoxes<dim>> & gdb_ext,
								 openfpm::vector<size_t> & gdb_ext_markers)
{
	size_t start = gdb_ext_markers.get(sub_id);
	size_t stop = gdb_ext_markers.get(sub_id+1);
	return get_gdb_ext(gdb_ext,start,stop,def_id);
}

/*! \brief It store the information about the local external ghost box
 *
 *
 */
template <unsigned int dim> struct e_lbox_id
{
	//! Box defining the external ghost box in local coordinates
	::Box<dim,long int> ebox;

	//! Has this external ghost box initialized
	bool initialized = false;

	//! sub-domain id of the non-extended sub-domain
	size_t sub;

	//! to which gdb_ext this external ghost box is linked with
	//! (in case of grid defined everywhere these two number match)
	size_t sub_gdb_ext;

	//! external ghost box linked to this internal ghost box
	size_t k;

	//! Sector position of the local external ghost box
	comb<dim> cmb;
};


/*! \brief Add a local external ghost box
 *
 * \param le_sub sub-domain id
 * \param s id of the external ghost box
 * \param j id of the domain not expanded
 * \param k id of the grid (gdb_ext) this external ghost box is linked with
 * \param bid where to add the local external ghost box
 * \param box the box to add
 * \param cmb quadrant in which the external ghost box live
 *
 */
template<unsigned int dim> inline void add_loc_eg_box(size_t le_sub,
											   size_t se,
											   size_t j,
											   size_t k,
		                                       openfpm::vector<e_lbox_id<dim>> & bid,
		                                       const Box<dim,long int> & ebox,
											   comb<dim> & cmb)
{
	bid.add();

	bid.last().ebox = ebox;

	bid.last().sub = se;
	bid.last().sub_gdb_ext = k;

	bid.last().cmb = cmb;
	bid.last().k = j;
	bid.last().initialized = true;
}

/*! Add an entry of for an external ghost box
 *
 * \param k sub-domain id to which this external ghost box is linked
 * \param cmb quadrant where the received linked internal ghost box live
 * \param output external ghost box
 * \param g_id global id of the external ghost box, in general this id is communicated
 *             by the processor that has the linked internal ghost-box
 * \param origin domain where this external ghost box is linked
 * \param p1 origin of the received internal ghost box
 *
 */
template<unsigned int dim> inline void add_eg_box(size_t k,
		                                   const comb<dim> & cmb,
										   const Box<dim,long int> & output,
										   size_t g_id,
										   const Point<dim,long int> & origin,
										   const Point<dim,long int> & p1,
										   openfpm::vector<e_box_id<dim>> & bid)
{
	// link

	size_t sub_id = k;

	e_box_id<dim> bid_t;
	bid_t.sub = sub_id;
	bid_t.cmb = cmb;
	bid_t.cmb.sign_flip();
	::Box<dim,long int> ib = output;
	bid_t.g_e_box = ib;
	bid_t.g_id = g_id;

	// Translate in local coordinate for gdb_ext
	Box<dim,long int> tb = ib;
	tb -= origin;
	bid_t.l_e_box = tb;

	// Translate in local coordinates for the received box
	Box<dim,long int> tbr = ib;
	bid_t.lr_e_box = tbr;

	bid.add(bid_t);
}

/*! \brief Result of the itersection of a box with an array of boxes
 *
 *
 */
template<unsigned int dim>
struct result_box
{
	//! id of the box in the array that produced an non-empty intersection
	size_t id;

	//! valid result of the itersection
	Box<dim,long int> bx;
};

/*! \brief Intersect a box with an array of boxes
 *
 * \param bx_def array of boxes
 * \param bx box to intersect with
 * \param use_bx_def in case of false the box bx is added to the result array
 *        and nothing is performed
 * \param result results of the intersections
 *
 */
template<unsigned int dim>
void bx_intersect(openfpm::vector<Box<dim,long int>> & bx_def,
		          bool use_bx_def,
				  Box<dim,long int> & bx,
				  openfpm::vector_std<result_box<dim>> & result)
{
	result.clear();

	if (use_bx_def == false)
	{
		result_box<dim> tmp;
		tmp.bx = bx;
		tmp.id = 0;

		result.add(tmp);
		return;
	}

	for (size_t i = 0 ; i < bx_def.size() ; i++)
	{
		result_box<dim> inte;
		if (bx.Intersect(bx_def.get(i),inte.bx))
		{
			inte.id = i;
			result.add(inte);
		}
	}
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

/*! \brief set of internal ghost box to send
 *
 *
 */
template<unsigned int dim>
struct ibox_send
{
	//! global id
	size_t gid;

	//! internal ghost box
	Box<dim,long int> ibox;
};



/*! \brief it store an internal ghost box, the linked external ghost box and the sub-domain from where
 *  it come from as internal ghost box
 *
 */
template<unsigned int dim> struct i_lbox_id
{
	//! Box
	::Box<dim,long int> box;

	//! sub-domain id (of the extended sub-domain).
	//! sub and sub_gdb_ext match if use_bx_def is not used in the
	//! distributed grid
	size_t sub;

	//! to which gdb_ext this internal ghost box is linked with
	size_t sub_gdb_ext;

	//! external ghost box linked to this internal ghost box
	openfpm::vector<size_t> k;

	//! combination
	comb<dim> cmb;
};


/*! \brief For each external ghost id, it contain a set of sub-domain at which this
 *         external box is linked
 *
 *
 */
template<unsigned int dim>
struct e_box_multi
{
	//! set sub-domain at which with external ghost is linked
	//! The eb_list are id for the eb_box list
	openfpm::vector<size_t> eb_list;

	//! This is the id in eb_list that contain an external ghost box
	//! able to store the full received box
	size_t full_match;

	//! It store the id of linked eg_box (eg_box.get(e_id) is valid, this mean that e_id store processor id information)
	//! \see ProctoID function
	size_t e_id;
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

	//! total number of received points
	size_t recv_pnt;

	//! Number of received boxes
	size_t n_r_box;
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
