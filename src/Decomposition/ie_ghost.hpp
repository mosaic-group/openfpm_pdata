/*
 * ie_ghost.hpp
 *
 *  Created on: Aug 8, 2015
 *      Author: i-bird
 */

#ifndef SRC_DECOMPOSITION_IE_GHOST_HPP_
#define SRC_DECOMPOSITION_IE_GHOST_HPP_

#include "common.hpp"

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
	void create_box_nn_processor_ext(Ghost<dim,T> & ghost, openfpm::vector<SpaceBox<dim,T>> & sub_domains)
	{
/*		box_nn_processor_int.resize(sub_domains.size());
		proc_int_box.resize(getNNProcessors());

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
				Box_dom & proc_int_box_g = proc_int_box.get(ProctoID(p_id));

				// get the set of sub-domains of the adjacent processor p_id
				openfpm::vector< ::Box<dim,T> > & nn_processor_subdomains_g = nn_processor_subdomains[p_id].bx;

				// near processor sub-domain intersections
				openfpm::vector< ::Box<dim,T> > & box_nn_processor_int_gg = box_nn_processor_int.get(i).get(j).bx;

				// for each near processor sub-domain intersect with the enlarged local sub-domain and store it
				for (size_t b = 0 ; b < nn_processor_subdomains_g.size() ; b++)
				{
					::Box<dim,T> bi;

					bool intersect = sub_with_ghost.Intersect(::Box<dim,T>(nn_processor_subdomains_g.get(b)),bi);

					if (intersect == true)
					{
						struct p_box pb;

						pb.box = bi;
						pb.proc = p_id;
						pb.lc_proc = ProctoID(p_id);

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
						size_t p_idp = ProctoID(p_id);
						for (k = 0 ; k < proc_adj_box.get(p_idp).size() ; k++)
						{
							if (proc_adj_box.get(p_idp).get(k) == i)
								break;
						}
						if (k == proc_adj_box.get(p_idp).size())
							std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " sub-domain not found\n";

						proc_int_box_g.ebx.last().id = (k * nn_processor_subdomains_g.size() + b) * v_cl.getProcessingUnits() + p_id;
					}
				}
			}
		}*/
	}

	/*! \brief Create the box_nn_processor_int (nbx part) structure, the geo_cell list and proc_int_box
	 *
	 * This structure store for each sub-domain of this processors the boxes that come from the intersection
	 * of the near processors sub-domains enlarged by the ghost size (Internal ghost box). These boxes
	 * fill a geometrical cell list. The proc_int_box store the same information ordered by near processors
	 *
	 * \param ghost margins
	 *
	 * \note Are the B8_0 B9_0 B9_1 B5_0 boxes in calculateGhostBoxes
	 * \see calculateGhostBoxes
	 *
	 */
	void create_box_nn_processor_int(Ghost<dim,T> & ghost, openfpm::vector<SpaceBox<dim,T>> & sub_domains, const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processors, const nn_prcs<dim,T> & nn_p )
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
				openfpm::vector< ::Box<dim,T> > & nn_p_box = nn_p.get nn_processor_subdomains[p_id].bx;

				// get the local processor id
				size_t lc_proc = nn_processor_subdomains[p_id].id;

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
					p_box b_int;
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
						Box_dom & pr_box_int = proc_int_box.get(ProctoID(p_id));
						Box_sub<dim,T> sb;
						sb = b_int.box;
						sb.sub = i;

						// Search for the correct id
						size_t s = 0;
						size_t p_idp = ProctoID(p_id);
						for (s = 0 ; s < proc_adj_box.get(p_idp).size() ; s++)
						{
							if (proc_adj_box.get(p_idp).get(s) == i)
								break;
						}
						if (s == proc_adj_box.get(p_idp).size())
							std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " sub-domain not found\n";

						sb.id = (k * proc_adj_box.get(p_idp).size() + s) * v_cl.getProcessingUnits() + v_cl.getProcessUnitID();

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
};


#endif /* SRC_DECOMPOSITION_IE_GHOST_HPP_ */
