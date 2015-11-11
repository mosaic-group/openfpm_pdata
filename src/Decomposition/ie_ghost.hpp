/*
 * ie_ghost.hpp
 *
 *  Created on: Aug 8, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_IE_GHOST_HPP_
#define SRC_DECOMPOSITION_IE_GHOST_HPP_

#include "common.hpp"
#include "nn_processor.hpp"

/*! \brief structure that store and compute the internal and external local ghost box
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 *
 * \see CartDecomposition
 *
 */
template<unsigned int dim, typename T>
class ie_ghost
{
	//! for each sub-domain (first vector), contain the list (nested vector) of the neighborhood processors
	//! and for each processor contain the boxes calculated from the intersection
	//! of the sub-domains + ghost with the near-by processor sub-domain () and the other way around
	//! \see calculateGhostBoxes
	openfpm::vector< openfpm::vector< Box_proc<dim,T> > > box_nn_processor_int;

	//! It store the same information of box_nn_processor_int organized by processor id
	openfpm::vector< Box_dom<dim,T> > proc_int_box;

	// External ghost boxes for this processor, indicated with G8_0 G9_0 ...
	openfpm::vector<p_box<dim,T> > vb_ext;

	// Internal ghost boxes for this processor domain, indicated with B8_0 B9_0 ..... in the figure
	// below as a linear vector
	openfpm::vector<p_box<dim,T> > vb_int;

	//! Cell-list that store the geometrical information of the internal ghost boxes
	CellList<dim,T,FAST> geo_cell;

protected:

	/*! \brief Initialize the geo cell list structure
	 *
	 * The geo cell list structure exist to speed up the labelling the points if they fall on some
	 * internal ghost
	 *
	 */
	void Initialize_geo_cell(const Box<dim,T> & domain, const size_t (&div)[dim] ,const Point<dim,T> & orig)
	{
		// Initialize the geo_cell structure
		geo_cell.Initialize(domain,div,orig);
	}

	/*! \brief Create the box_nn_processor_int (bx part)  structure
	 *
	 * This structure store for each sub-domain of this processors enlarged by the ghost size the boxes that
	 *  come from the intersection with the near processors sub-domains (External ghost box)
	 *
	 * \param ghost margins
	 *
	 * \note Are the G8_0 G9_0 G9_1 G5_0 boxes in calculateGhostBoxes
	 * \see calculateGhostBoxes
	 *
	 */
	void create_box_nn_processor_ext(Vcluster & v_cl, Ghost<dim,T> & ghost, openfpm::vector<SpaceBox<dim,T>> & sub_domains, const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processor, const nn_prcs<dim,T> & nn_p)
	{
		box_nn_processor_int.resize(sub_domains.size());
		proc_int_box.resize(nn_p.getNNProcessors());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			SpaceBox<dim,T> sub_with_ghost = sub_domains.get(i);

			// enlarge the sub-domain with the ghost
			sub_with_ghost.enlarge(ghost);

			// resize based on the number of adjacent processors
			box_nn_processor_int.get(i).resize(box_nn_processor.get(i).size());

			// For each processor adjacent to this sub-domain
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// Contiguous processor
				size_t p_id = box_nn_processor.get(i).get(j);

				// store the box in proc_int_box storing from which sub-domain they come from
				Box_dom<dim,T> & proc_int_box_g = proc_int_box.get(nn_p.ProctoID(p_id));

				// get the set of sub-domains of the adjacent processor p_id
				const openfpm::vector< ::Box<dim,T> > & nn_processor_subdomains_g = nn_p.getExternalAdjSubdomain(p_id).bx;

				// near processor sub-domain intersections
				openfpm::vector< ::Box<dim,T> > & box_nn_processor_int_gg = box_nn_processor_int.get(i).get(j).bx;

				// for each near processor sub-domain intersect with the enlarged local sub-domain and store it
				for (size_t b = 0 ; b < nn_processor_subdomains_g.size() ; b++)
				{
					::Box<dim,T> bi;

					bool intersect = sub_with_ghost.Intersect(::Box<dim,T>(nn_processor_subdomains_g.get(b)),bi);

					if (intersect == true)
					{
						struct p_box<dim,T> pb;

						pb.box = bi;
						pb.proc = p_id;
						pb.lc_proc = nn_p.ProctoID(p_id);

						//
						// Updating
						//
						// vb_ext
						// box_nn_processor_int
						// proc_int_box
						//
						// They all store the same information but organized in different ways
						// read the description of each for more information
						//
						vb_ext.add(pb);
						box_nn_processor_int_gg.add(bi);
						proc_int_box_g.ebx.add();
						proc_int_box_g.ebx.last() = bi;
						proc_int_box_g.ebx.last().sub = i;

						// Search for the correct id
						size_t k = 0;
						size_t p_idp = nn_p.ProctoID(p_id);
						for (k = 0 ; k < nn_p.getInternalAdjSubdomain(p_idp).size() ; k++)
						{
							if (nn_p.getInternalAdjSubdomain(p_idp).get(k) == i)
								break;
						}
						if (k == nn_p.getInternalAdjSubdomain(p_idp).size())
							std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " sub-domain not found\n";

						proc_int_box_g.ebx.last().id = (k * nn_processor_subdomains_g.size() + b) * v_cl.getProcessingUnits() + p_id;
					}
				}
			}
		}
	}

	/*! \brief Create the box_nn_processor_int (nbx part) structure, the geo_cell list and proc_int_box
	 *
	 * This structure store for each sub-domain of this processors the boxes that come from the intersection
	 * of the near processors sub-domains enlarged by the ghost size (Internal ghost box). These boxes
	 * fill a geometrical cell list. The proc_int_box store the same information ordered by near processors
	 *
	 * \param v_cl Virtual cluster
	 * \param ghost margins
	 * \param sub_domains
	 * \param box_nn_processors sub-domains of the adjacent processors
	 * \param nn_prcs structure that store the adjacent processor information
	 * \param geo_cell Cell list that store the subdomain information
	 *
	 * \note Are the B8_0 B9_0 B9_1 B5_0 boxes in calculateGhostBoxes
	 * \see calculateGhostBoxes
	 *
	 */
	void create_box_nn_processor_int(Vcluster & v_cl, Ghost<dim,T> & ghost, openfpm::vector<SpaceBox<dim,T>> & sub_domains, const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processor, const nn_prcs<dim,T> & nn_p)
	{
		box_nn_processor_int.resize(sub_domains.size());
		proc_int_box.resize(nn_p.getNNProcessors());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			// For each processor contiguous to this sub-domain
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// Contiguous processor
				size_t p_id = box_nn_processor.get(i).get(j);

				// get the set of sub-domains of the contiguous processor p_id
				const openfpm::vector< ::Box<dim,T> > & nn_p_box = nn_p.getExternalAdjSubdomain(p_id).bx;

				// get the local processor id
				size_t lc_proc = nn_p.getAdjacentProcessor(p_id);

				// For each near processor sub-domains enlarge and intersect with the local sub-domain and store the result
				for (size_t k = 0 ; k < nn_p_box.size() ; k++)
				{

					// enlarge the near-processor sub-domain
					::Box<dim,T> n_sub = nn_p_box.get(k);

					// local sub-domain
					::SpaceBox<dim,T> l_sub = sub_domains.get(i);

					// Create a margin of ghost size around the near processor sub-domain
					n_sub.enlarge(ghost);

					// Intersect with the local sub-domain
					p_box<dim,T> b_int;
					bool intersect = n_sub.Intersect(l_sub,b_int.box);

					// store if it intersect
					if (intersect == true)
					{
						// the box fill with the processor id
						b_int.proc = p_id;

						// fill the local processor id
						b_int.lc_proc = lc_proc;

						//
						// Updating
						//
						// vb_int
						// box_nn_processor_int
						// proc_int_box
						//
						// They all store the same information but organized in different ways
						// read the description of each for more information
						//

						// add the box to the near processor sub-domain intersections
						openfpm::vector< ::Box<dim,T> > & p_box_int = box_nn_processor_int.get(i).get(j).nbx;
						p_box_int.add(b_int.box);
						vb_int.add(b_int);

						// store the box in proc_int_box storing from which sub-domain they come from
						Box_dom<dim,T> & pr_box_int = proc_int_box.get(nn_p.ProctoID(p_id));
						Box_sub<dim,T> sb;
						sb = b_int.box;
						sb.sub = i;

						// Search for the correct id
						size_t s = 0;
						size_t p_idp = nn_p.ProctoID(p_id);
						for (s = 0 ; s < nn_p.getInternalAdjSubdomain(p_idp).size() ; s++)
						{
							if (nn_p.getInternalAdjSubdomain(p_idp).get(s) == i)
								break;
						}
						if (s == nn_p.getInternalAdjSubdomain(p_idp).size())
							std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " sub-domain not found\n";

						sb.id = (k * nn_p.getInternalAdjSubdomain(p_idp).size() + s) * v_cl.getProcessingUnits() + v_cl.getProcessUnitID();

						pr_box_int.ibx.add(sb);

						// update the geo_cell list

						// get the cells this box span
						const grid_key_dx<dim> p1 = geo_cell.getCellGrid(b_int.box.getP1());
						const grid_key_dx<dim> p2 = geo_cell.getCellGrid(b_int.box.getP2());

						// Get the grid and the sub-iterator
						auto & gi = geo_cell.getGrid();
						grid_key_dx_iterator_sub<dim> g_sub(gi,p1,p2);

						// add the box-id to the cell list
						while (g_sub.isNext())
						{
							auto key = g_sub.get();
							geo_cell.addCell(gi.LinId(key),vb_int.size()-1);
							++g_sub;
						}
					}
				}
			}
		}
	}

public:

	/*! \brief Get the number of Internal ghost boxes for one processor
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \return the number of internal ghost
	 *
	 */
	inline size_t getProcessorNIGhost(size_t id) const
	{
		return proc_int_box.get(id).ibx.size();
	}

	/*! \brief Get the number of External ghost boxes for one processor id
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \return the number of external ghost
	 *
	 */
	inline size_t getProcessorNEGhost(size_t id) const
	{
		return proc_int_box.get(id).ebx.size();
	}

	/*! \brief Get the j Internal ghost box for one processor
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one internal ghost box)
	 * \return the box
	 *
	 */
	inline const ::Box<dim,T> & getProcessorIGhostBox(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ibx.get(j);
	}

	/*! \brief Get the j External ghost box
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one external ghost box)
	 * \return the box
	 *
	 */
	inline const ::Box<dim,T> & getProcessorEGhostBox(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ebx.get(j);
	}

	/*! \brief Get the j Internal ghost box id
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one internal ghost box)
	 * \return the box
	 *
	 */
	inline size_t getProcessorIGhostId(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ibx.get(j).id;
	}

	/*! \brief Get the j External ghost box id
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one external ghost box)
	 * \return the box
	 *
	 */
	inline size_t getProcessorEGhostId(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ebx.get(j).id;
	}


	/*! \brief Get the local sub-domain at witch belong the internal ghost box
	 *
	 * \param id adjacent processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one internal ghost box)
	 * \return sub-domain at which belong the internal ghost box
	 *
	 */
	inline size_t getProcessorIGhostSub(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ibx.get(j).sub;
	}

	/*! \brief Get the local sub-domain at witch belong the external ghost box
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one external ghost box)
	 * \return sub-domain at which belong the external ghost box
	 *
	 */
	inline size_t getProcessorEGhostSub(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ebx.get(j).sub;
	}

	/*! \brief Return the total number of the calculated internal ghost boxes
	 *
	 * \return the number of internal ghost boxes
	 *
	 */
	inline size_t getNIGhostBox() const
	{
		return vb_int.size();
	}

	/*! \brief Given the internal ghost box id, it return the internal ghost box
	 *
	 * \return the internal ghost box
	 *
	 */
	inline ::Box<dim,T> getIGhostBox(size_t b_id) const
	{
		return vb_int.get(b_id).box;
	}

	/*! \brief Given the internal ghost box id, it return the near processor at witch belong
	 *         or the near processor that produced this internal ghost box
	 *
	 * \return the processor id of the ghost box
	 *
	 */
	inline size_t getIGhostBoxProcessor(size_t b_id) const
	{
		return vb_int.get(b_id).proc;
	}

	/*! \brief Get the number of the calculated external ghost boxes
	 *
	 * \return the number of external ghost boxes
	 *
	 */
	inline size_t getNEGhostBox() const
	{
		return vb_ext.size();
	}

	/*! \brief Given the external ghost box id, it return the external ghost box
	 *
	 * \return the external ghost box
	 *
	 */
	inline ::Box<dim,T> getEGhostBox(size_t b_id) const
	{
		return vb_ext.get(b_id).box;
	}

	/*! \brief Given the external ghost box id, it return the near processor at witch belong
	 *         or the near processor that produced this external ghost box
	 *
	 * \return the processor id of the external ghost box
	 *
	 */
	inline size_t getEGhostBoxProcessor(size_t b_id) const
	{
		return vb_ext.get(b_id).proc;
	}

	/*! /brief Given a point it return the set of boxes in which the point fall
	 *
	 * \param p Point to check
	 * \return An iterator with the id's of the internal boxes in which the point fall
	 *
	 */
	auto getInternalIDBoxes(Point<dim,T> & p) -> decltype(geo_cell.getIterator(geo_cell.getCell(p)))
	{
		return geo_cell.getIterator(geo_cell.getCell(p));
	}

	/*! \brief if the point fall into the ghost of some near processor it return the processors id's in which
	 *  it fall
	 *
	 * \param p Point
	 * \return iterator of the processors id's
	 *
	 */
	inline auto labelPoint(Point<dim,T> & p) -> decltype(geo_cell.getIterator(geo_cell.getCell(p)))
	{
		return geo_cell.getIterator(geo_cell.getCell(p));
	}

	openfpm::vector<size_t> ids;

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 * (Internal ghost)
	 *
	 * \tparam id type of if to get box_id processor_id lc_processor_id
	 * \param p Particle position
	 * \param opt intersection boxes of the same processor can overlap, so in general the function
	 *        can produce more entry with the same processor, the UNIQUE option eliminate double entries
	 *        (UNIQUE) is for particle data (MULTIPLE) is for grid data [default MULTIPLE]
	 *
	 * \param return the processor ids
	 *
	 */
	template <typename id> inline const openfpm::vector<size_t> ghost_processorID(Point<dim,T> & p, const int opt = MULTIPLE)
	{
		ids.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getIterator(geo_cell.getCell(p));

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInside(p) == true)
			{
				ids.add(id::id(vb_int.get(bid),bid));
			}

			++cell_it;
		}

		// Make the id unique
		if (opt == UNIQUE)
			ids.unique();

		return ids;
	}

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 * (Internal ghost)
	 *
	 * \tparam id type of if to get box_id processor_id lc_processor_id
	 * \param p Particle position
	 *
	 * \param return the processor ids
	 *
	 */
	template<typename id, typename Mem> inline const openfpm::vector<size_t> ghost_processorID(const encapc<1,Point<dim,T>,Mem> & p, const int opt = MULTIPLE)
	{
		ids.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getIterator(geo_cell.getCell(p));

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInside(p) == true)
			{
				ids.add(id::id(vb_int.get(bid),bid));
			}

			++cell_it;
		}

		// Make the id unique
		if (opt == UNIQUE)
			ids.unique();

		return ids;
	}

	/*! \brief write the information about the ghost in vtk format
	 *
	 * 1) internal_ghost_X.vtk Internal ghost boxes for the local processor (X)
	 * 2) external_ghost_X.vtk External ghost boxes for the local processor (X)
	 *
	 * \param output directory
	 * \param p_id processor rank
	 *
	 */
	bool write(std::string output, size_t p_id) const
	{
		//! internal_ghost_X.vtk Internal ghost boxes for the local processor (X)
		VTKWriter<openfpm::vector<::Box<dim,T>>,VECTOR_BOX> vtk_box3;
		for (size_t p = 0 ; p < box_nn_processor_int.size() ; p++)
		{
			for (size_t s = 0 ; s < box_nn_processor_int.get(p).size() ; s++)
			{
				vtk_box3.add(box_nn_processor_int.get(p).get(s).nbx);
			}
		}
		vtk_box3.write(output + std::string("internal_ghost_") + std::to_string(p_id) + std::string(".vtk"));

		//! external_ghost_X.vtk External ghost boxes for the local processor (X)
		VTKWriter<openfpm::vector<::Box<dim,T>>,VECTOR_BOX> vtk_box4;
		for (size_t p = 0 ; p < box_nn_processor_int.size() ; p++)
		{
			for (size_t s = 0 ; s < box_nn_processor_int.get(p).size() ; s++)
			{
				vtk_box4.add(box_nn_processor_int.get(p).get(s).bx);
			}
		}
		vtk_box4.write(output + std::string("external_ghost_") + std::to_string(p_id) + std::string(".vtk"));

		return true;
	}

	/*! \brief Check if the ie_loc_ghosts contain the same information
	 *
	 * \param ele Element to check
	 *
	 */
	bool is_equal(ie_ghost<dim,T> & ig)
	{
		if (getNEGhostBox() != ig.getNEGhostBox())
			return false;

		if (getNIGhostBox() != ig.getNIGhostBox())
			return false;

		for (size_t i = 0 ; i < getNIGhostBox() ; i++)
		{
			if (getProcessorNIGhost(i) != ig.getProcessorNIGhost(i))
				return false;
			for (size_t j = 0 ; j < getProcessorNIGhost(i) ; j++)
			{
				if (getProcessorIGhostBox(i,j) != ig.getProcessorIGhostBox(i,j))
					return false;
				if (getProcessorIGhostId(i,j) != ig.getProcessorIGhostId(i,j))
					return false;
				if (getProcessorIGhostSub(i,j) != ig.getProcessorIGhostSub(i,j))
					return false;
			}
			if (getIGhostBox(i) != ig.getIGhostBox(i))
				return false;
			if (getIGhostBoxProcessor(i) != ig.getIGhostBoxProcessor(i))
				return false;
		}

		for (size_t i = 0 ; i < getNEGhostBox() ; i++)
		{
			if (getProcessorNEGhost(i) != ig.getProcessorNEGhost(i))
				return false;
			for (size_t j = 0 ; j < getProcessorNEGhost(i) ; j++)
			{
				if (getProcessorEGhostBox(i,j) != ig.getProcessorEGhostBox(i,j))
					return false;
				if (getProcessorEGhostId(i,j) != ig.getProcessorEGhostId(i,j))
					return false;
				if (getProcessorEGhostSub(i,j) != ig.getProcessorEGhostSub(i,j))
					return false;
			}
			if (getEGhostBox(i) != ig.getEGhostBox(i))
				return false;
			if (getEGhostBoxProcessor(i) != ig.getEGhostBoxProcessor(i))
				return false;
		}
	}

	/*! \brief Check if the ie_loc_ghosts contain the same information with the exception of the ghost part
	 * It is anyway required that the ghost come from the same sub-domains decomposition
	 *
	 * \param ele Element to check
	 *
	 */
	bool is_equal_ng(ie_ghost<dim,T> & ig)
	{
		Box<dim,T> bt;

		if (getNEGhostBox() != ig.getNEGhostBox())
			return false;

		if (getNIGhostBox() != ig.getNIGhostBox())
			return false;

		for (size_t i = 0 ; i < getNIGhostBox() ; i++)
		{
			if (getProcessorNIGhost(i) != ig.getProcessorNIGhost(i))
				return false;
			for (size_t j = 0 ; j < getProcessorNIGhost(i) ; j++)
			{
				if (getProcessorIGhostBox(i,j).intersect(ig.getProcessorIGhostBox(i,j),bt) == false)
					return false;
				if (getProcessorIGhostId(i,j).intersect(ig.getProcessorIGhostId(i,j),bt) == false)
					return false;
				if (getProcessorIGhostSub(i,j) != ig.getProcessorIGhostSub(i,j))
					return false;
			}
			if (getIGhostBox(i) != ig.getIGhostBox(i))
				return false;
			if (getIGhostBoxProcessor(i) != ig.getIGhostBoxProcessor(i))
				return false;
		}

		for (size_t i = 0 ; i < getNEGhostBox() ; i++)
		{
			if (getProcessorNEGhost(i) != ig.getProcessorNEGhost(i))
				return false;
			for (size_t j = 0 ; j < getProcessorNEGhost(i) ; j++)
			{
				if (getProcessorEGhostBox(i,j).intersect(ig.getProcessorEGhostBox(i,j),bt) == false)
					return false;
				if (getProcessorEGhostId(i,j),intersect(ig.getProcessorEGhostId(i,j),bt) == false)
					return false;
				if (getProcessorEGhostSub(i,j) != ig.getProcessorEGhostSub(i,j))
					return false;
			}
			if (getEGhostBox(i) != ig.getEGhostBox(i))
				return false;
			if (getEGhostBoxProcessor(i).intersect(ig.getEGhostBoxProcessor(i),bt) == false)
				return false;
		}
	}
};


#endif /* SRC_DECOMPOSITION_IE_GHOST_HPP_ */
