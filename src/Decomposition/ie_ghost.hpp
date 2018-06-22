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
#include "Decomposition/shift_vect_converter.hpp"


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

	//! External ghost boxes for this processor
	openfpm::vector<p_box<dim,T> > vb_ext;

	//! Internal ghost boxes for this processor domain
	openfpm::vector<p_box<dim,T> > vb_int;

	//! Cell-list that store the geometrical information of the internal ghost boxes
	CellList<dim,T,Mem_fast<>,shift<dim,T>> geo_cell;

	//! shift vectors
	openfpm::vector<Point<dim,T>> shifts;

	//! Temporal buffers to return temporal information for ghost_processorID
	openfpm::vector<std::pair<size_t,size_t>> ids_p;

	//! Temporal buffers to return temporal information
	openfpm::vector<size_t> ids;

	//! shift converter
	shift_vect_converter<dim,T> sc_convert;

	/*! \brief Given a local sub-domain i, it give the id of such sub-domain in the sent list
	 *         for the processor p_id
	 *
	 * Processor 5 send its sub-domains to processor 6 and will receive the list from 6
	 *
	 * This function search if a local sub-domain has been sent to a processor p_id, if
	 * found it return at witch position is in the list of the sent sub-domains
	 *
	 * \param nn_p structure that store the processor graph as near processor
	 * \param p_id near processor rank
	 * \param i sub-domain
	 *
	 * \return Given a local sub-domain i, it give the id of such sub-domain in the sent list
	 *         for the processor p_id
	 *
	 */
	inline size_t link_ebx_ibx(const nn_prcs<dim,T> & nn_p, size_t p_id, size_t i)
	{
		// Search for the correct id
		size_t k = 0;
		size_t p_idp = nn_p.ProctoID(p_id);
		for (k = 0 ; k < nn_p.getSentSubdomains(p_idp).size() ; k++)
		{
			if (nn_p.getSentSubdomains(p_idp).get(k) == i)
				break;
		}
		if (k == nn_p.getSentSubdomains(p_idp).size())
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " sub-domain not found\n";

		return k;
	}

	/*! \brief This is the external and internal ghost box link formula
	*
	* This formula is pretty important and require an extensive explanation
	*
	* \verbatim

	+------------+
	|            |
	|            +---+---------+
	| Processor 5|   |         |
	|            | E | Proc 6  |
	|  Sub 0     | 0 |         |
	|            | _ | Sub 9   |
	|            | 9 |         |
	|            |   |         |
	|            +---+---------+
	|            |
	+------------+

	* \endverbatim
	*
	* E0_6 is an external ghost box from the prospective of processor 5 and an internal
	* ghost boxes from the prospective of processor 6. So for every external
	* ghost box that processor 5 compute, exist an internal ghost box in processor 6
	*
	* Here we link this information with an unique id, for processor 5 and 6.
	* Consider Processor 5 sending to processor 6
	* its sub-domains, including the one in figure with id 0 in the list, and
	* receive from processor 6 the sub-domain in figure as id 9. Consider also
	*  we have 16 processor. E0_9 come from the intersection of the expanded sub-domain
	* 0 with 9 (Careful the id is related to the send and receive position in the list)
	* and the intersection is in the sector 0
	*
	*
	* The id of the external box (for processor 5) is calculated as
	*
	* ((k * N_b + b) * v_cl.getProcessingUnits() + p_id) * openfpm::math::pow(3,dim) + c.lin()
	*
	* The parameter assume a different meaning if they the formula is used for calculating
	* external/internal ghost boxes id
	*
	* \param k expanded sub-domain sent/received to/from p_id ( 0 )
	* \param b sub-domain received/sent from/to p_id ( 9 )
	* \param p_id processor id ( 6 )
	* \param c sector where the sub-domain b live
	* \param N_b number of sub-domain received/sent from/to p_id
	* \param v_cl Vcluster
	* \param ei indicate if the formula is used to calculate external (true) or internal (false) ids
	*
	* \return id of the external/internal ghost
	*
	* \note To an explanation about the sectors see getShiftVectors
	*
	*/
	inline size_t ebx_ibx_form(size_t k, size_t b, size_t p_id, const comb<dim> & c ,size_t N_b, Vcluster & v_cl, const bool ei)
	{
		comb<dim> cext = c;

		if (ei == true)
			cext.sign_flip();

		return ((k * N_b + b) * v_cl.getProcessingUnits() + p_id) * openfpm::math::pow(3,dim) + cext.lin();
	}

protected:

	/*! \brief Here we generare the shift vectors
	 *
	 * \param domain box that describe the domain
	 *
	 */
	void generateShiftVectors(const Box<dim,T> & domain, size_t (& bc)[dim])
	{
		sc_convert.generateShiftVectors(domain,bc,shifts);
	}

	/*! \brief Initialize the geo cell list structure
	 *
	 * The geo cell list structure exist to speed up the labelling the points if they fall on some
	 * internal ghost
	 *
	 * \param domain where the cell list is defined
	 * \param div number of division of the cell list
	 *
	 */
	void Initialize_geo_cell(const Box<dim,T> & domain, const size_t (&div)[dim])
	{
		// Initialize the geo_cell structure
		geo_cell.Initialize(domain,div,0);
	}

	/*! \brief Create the box_nn_processor_int (bx part)  structure
	 *
	 * For each sub-domain of the local processor it store the intersection between the enlarged
	 * sub-domain of the calling processor with the adjacent processors sub-domains (External ghost box)
	 *
	 * \param v_cl Virtual cluster
	 * \param ghost margins
	 * \param sub_domains vector of local sub-domains
	 * \param box_nn_processor it will store for each sub-domain the near processors
	 * \param nn_p contain the sub-domains of the near processors
	 *
	 * \note Are the G8_0 G9_0 G9_1 G5_0 boxes in calculateGhostBoxes
	 * \see calculateGhostBoxes
	 *
	 */
	void create_box_nn_processor_ext(Vcluster & v_cl,
			                         Ghost<dim,T> & ghost,
									 openfpm::vector<SpaceBox<dim,T>> & sub_domains,
									 const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processor,
									 const nn_prcs<dim,T> & nn_p)
	{
		box_nn_processor_int.resize(sub_domains.size());
		proc_int_box.resize(nn_p.getNNProcessors());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			SpaceBox<dim,T> sub_with_ghost = sub_domains.get(i);

			// enlarge the sub-domain with the ghost
			sub_with_ghost.enlarge(ghost);

			// resize based on the number of near processors
			box_nn_processor_int.get(i).resize(box_nn_processor.get(i).size());

			// For each processor near to this sub-domain
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// near processor
				size_t p_id = box_nn_processor.get(i).get(j);

				// used later
				Box_dom<dim,T> & proc_int_box_g = proc_int_box.get(nn_p.ProctoID(p_id));

				// Number of received sub-domains
				size_t n_r_sub = nn_p.getNRealSubdomains(p_id);

				// get the set of sub-domains, sector position, and real sub-domain id of the near processor p_id
				const openfpm::vector< ::Box<dim,T> > & nn_processor_subdomains_g = nn_p.getNearSubdomains(p_id);
				const openfpm::vector< comb<dim> > & nnpsg_pos = nn_p.getNearSubdomainsPos(p_id);
				const openfpm::vector< size_t > & r_sub = nn_p.getNearSubdomainsRealId(p_id);

				// used later
				openfpm::vector< ::Box<dim,T> > & box_nn_processor_int_gg = box_nn_processor_int.get(i).get(j).bx;

				// for each near processor sub-domain intersect with the enlarged local sub-domain and store it
				for (size_t b = 0 ; b < nn_processor_subdomains_g.size() ; b++)
				{
					::Box<dim,T> bi;
					::Box<dim,T> sub_bb(nn_processor_subdomains_g.get(b));

					bool intersect = sub_with_ghost.Intersect(sub_bb,bi);

					if (intersect == true)
					{
						struct p_box<dim,T> pb;

						pb.box = bi;
						pb.proc = p_id;
						pb.lc_proc = nn_p.ProctoID(p_id);
						pb.shift_id = (size_t)-1;

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
						proc_int_box_g.ebx.last().bx = bi;
						proc_int_box_g.ebx.last().sub = i;
						proc_int_box_g.ebx.last().cmb = nnpsg_pos.get(b);

						// Search where the sub-domain i is in the sent list for processor p_id
						size_t k = link_ebx_ibx(nn_p,p_id,i);

						proc_int_box_g.ebx.last().id = ebx_ibx_form(k,r_sub.get(b),p_id,nnpsg_pos.get(b),n_r_sub,v_cl,true);
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
	 * \param box_nn_processor sub-domains of the near processors
	 * \param nn_p structure that store the near processor sub-domains
	 *
	 * \note Are the B8_0 B9_0 B9_1 B5_0 boxes in calculateGhostBoxes
	 * \see calculateGhostBoxes
	 *
	 */
	void create_box_nn_processor_int(Vcluster & v_cl,
			                         Ghost<dim,T> & ghost,
									 openfpm::vector<SpaceBox<dim,T>> & sub_domains,
									 const openfpm::vector<openfpm::vector<long unsigned int> > & box_nn_processor,
									 const nn_prcs<dim,T> & nn_p)
	{
		box_nn_processor_int.resize(sub_domains.size());
		proc_int_box.resize(nn_p.getNNProcessors());

		// For each sub-domain
		for (size_t i = 0 ; i < sub_domains.size() ; i++)
		{
			// For each processor contiguous to this sub-domain
			for (size_t j = 0 ; j < box_nn_processor.get(i).size() ; j++)
			{
				// Near processor
				size_t p_id = box_nn_processor.get(i).get(j);

				// get the set of sub-domains of the near processor p_id
				const openfpm::vector< ::Box<dim,T> > & nn_p_box = nn_p.getNearSubdomains(p_id);

				// get the sector position for each sub-domain in the list
				const openfpm::vector< comb<dim> > nn_p_box_pos = nn_p.getNearSubdomainsPos(p_id);

				// get the real sub-domain id for each sub-domain
				const openfpm::vector<size_t> r_sub = nn_p.getNearSubdomainsRealId(p_id);

				// get the local processor id
				size_t lc_proc = nn_p.getNearProcessor(p_id);

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

						// fill the shift id
						b_int.shift_id = convertShift(nn_p_box_pos.get(k));

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
						Box_sub<dim,T> sb;
						sb.bx = b_int.box;
						sb.sub = i;
						sb.r_sub = r_sub.get(k);
						sb.cmb = nn_p_box_pos.get(k);

						size_t p_idp = nn_p.ProctoID(p_id);

						// Search where the sub-domain i is in the sent list for processor p_id
						size_t s = link_ebx_ibx(nn_p,p_id,i);

						// calculate the id of the internal box
						sb.id = ebx_ibx_form(r_sub.get(k),s,v_cl.getProcessUnitID(),nn_p_box_pos.get(k),nn_p.getSentSubdomains(p_idp).size(),v_cl,false);

						Box_dom<dim,T> & pr_box_int = proc_int_box.get(nn_p.ProctoID(p_id));
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

	//! Default constructor
	ie_ghost() {};

	//! Copy constructor
	ie_ghost(const ie_ghost<dim,T> & ie)
	{
		this->operator =(ie);
	}

	//! Copy constructor
	ie_ghost(ie_ghost<dim,T> && ie)
	{
		this->operator=(ie);
	}

	//! Copy operator
	inline ie_ghost<dim,T> & operator=(ie_ghost<dim,T> && ie)
	{
		box_nn_processor_int.swap(ie.box_nn_processor_int);
		proc_int_box.swap(ie.proc_int_box);
		vb_ext.swap(ie.vb_ext);
		vb_int.swap(ie.vb_int);
		geo_cell.swap(ie.geo_cell);
		shifts.swap(ie.shifts);
		ids_p.swap(ie.ids_p);
		ids.swap(ie.ids);

		return *this;
	}

	//! Copy operator
	inline ie_ghost<dim,T> & operator=(const ie_ghost<dim,T> & ie)
	{
		box_nn_processor_int = ie.box_nn_processor_int;
		proc_int_box = ie.proc_int_box;
		vb_ext = ie.vb_ext;
		vb_int = ie.vb_int;
		geo_cell = ie.geo_cell;
		shifts = ie.shifts;
		ids_p = ie.ids_p;
		ids = ie.ids;

		return *this;
	}

	/*! It return the shift vector
	 *
	 * Consider a domain with some ghost, at the border of the domain the
	 * ghost must be treated in a special way, depending on the periodicity
	 * of the boundary
	 *
		\verbatim

															[1,1]
			+---------+------------------------+---------+
			| (1,-1)  |                        | (1,1)   |
			|   |     |    (1,0) --> 7         |   |     |
			|   v     |                        |   v     |
			|   6     |                        |   8     |
			+--------------------------------------------+
			|         |                        |         |
			|         |                        |         |
			|         |                        |         |
			| (-1,0)  |                        | (1,0)   |
			|    |    |                        |   |     |
			|    v    |      (0,0) --> 4       |   v     |
			|    3    |                        |   5     |
			|         |                        |         |
		 B	|         |                        |     A   |
		*	|         |                        |    *    |
			|         |                        |         |
			|         |                        |         |
			|         |                        |         |
			+--------------------------------------------+
			| (-1,-1) |                        | (-1,1)  |
			|    |    |   (-1,0) --> 1         |    |    |
			|    v    |                        |    v    |
			|    0    |                        |    2    |
			+---------+------------------------+---------+


		\endverbatim
	 *
	 *
	 * if a particle is bound in (1,0) linearized to 5, before communicate this particle (A in figure)
	 * must be shifted on -1.0 on x (B in figure)
	 *
	 * This function return the set of shift vectors that determine such shift, for example
	 * in the example above the shift at position 5 will be (0,-1.0)
	 *
	 * \return the shift vectors
	 *
	 */
	const openfpm::vector<Point<dim,T>> & getShiftVectors()
	{
		return shifts;
	}

	/*! It return the converted shift vector
	 *
	 * In high dimensions the number of shifts vectors explode exponentially, so we are
	 * expecting that some of the boundary is non periodic to reduce the numbers of shift
	 * vectors
	 *
	 * \return the shift vectors
	 *
	 */
	size_t convertShift(const comb<dim> & cmb)
	{
		return sc_convert.linId(cmb);
	}

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
		return proc_int_box.get(id).ibx.get(j).bx;
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
		return proc_int_box.get(id).ebx.get(j).bx;
	}

	/*! \brief Get the j External ghost box sector
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one external ghost box)
	 * \return the sector
	 *
	 */
	inline const comb<dim> & getProcessorEGhostPos(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ebx.get(j).cmb;
	}

	/*! \brief Get the ghost box sector of the external ghost box linked with the j internal ghost box
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one internal ghost box)
	 * \return the sector
	 *
	 */
	inline const comb<dim> & getProcessorIGhostPos(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ibx.get(j).cmb;
	}

	/*! \brief Get the j Internal ghost box id
	 *
	 * Every internal ghost box has a linked external ghost box, because they overlap
	 * and they must contain the same information (Think on a ghost_get). So if exist
	 *  an internal ghost box with id x, exist also an external ghost box with id x
	 *
	 * \param id near processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one internal ghost box)
	 * \return the box id
	 *
	 */
	inline size_t getProcessorIGhostId(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ibx.get(j).id;
	}

	/*! \brief Get the j External ghost box id
	 *
	 * Every external ghost box has a linked internal ghost box, because they overlap
	 * and they must contain the same information (Think on a ghost_get). So if exist
	 *  an internal ghost box with id x, exist also an external ghost box with id x
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

	/*! \brief Get the sub-domain send-id at witch belong the internal ghost box
	 *
	 * The internal ghost box is create from the intersection a local sub-domain
	 * and an extended sub-domain communicated from another processor. This function
	 * return the id of the sub-domain in the receiving list
	 *
	 * \param id adjacent processor list id (the id go from 0 to getNNProcessor())
	 * \param j box (each near processor can produce more than one internal ghost box)
	 * \return sub-domain at which belong the internal ghost box
	 *
	 */
	inline size_t getProcessorIGhostSSub(size_t id, size_t j) const
	{
		return proc_int_box.get(id).ibx.get(j).r_sub;
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
	 * \param b_id internal ghost box id
	 *
	 * \return the internal ghost box
	 *
	 */
	inline const ::Box<dim,T> & getIGhostBox(size_t b_id) const
	{
		return vb_int.get(b_id).box;
	}

	/*! \brief Given the internal ghost box id, it return the near processor at witch belong
	 *         or the near processor that produced this internal ghost box
	 *
	 * \param b_id internal ghost box id
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
	 * \param b_id external ghost box id
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
	 * \param b_id external ghost box id
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
	 *
	 * \return An iterator with the id's of the internal boxes in which the point fall
	 *
	 */
	auto getInternalIDBoxes(Point<dim,T> & p) -> decltype(geo_cell.getCellIterator(geo_cell.getCell(p)))
	{
		return geo_cell.getCellIterator(geo_cell.getCell(p));
	}

	/*! \brief if the point fall into the ghost of some near processor it return the processors id's in which
	 *  it fall
	 *
	 * \param p Point
	 * \return iterator of the processors id's
	 *
	 */
	inline auto labelPoint(Point<dim,T> & p) -> decltype(geo_cell.getCellIterator(geo_cell.getCell(p)))
	{
		return geo_cell.getCellIterator(geo_cell.getCell(p));
	}

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 * (Internal ghost)
	 *
	 * if the particle come from an internal ghost from the periodicity of the domain, position must be shifted
	 * this function return the id of the shift vector
	 *
	 * \see getShiftVector
	 *
	 * \tparam id type of id to get box_id processor_id lc_processor_id shift_id
	 *
	 * \param p Particle position
	 * \param opt intersection boxes of the same processor can overlap, so in general the function
	 *        can produce more entry with the same processor, the UNIQUE option eliminate double entries
	 *        (UNIQUE) is for particle data (MULTIPLE) is for grid data [default MULTIPLE]
	 * \return return the processor ids (not the rank, the id in the near processor list)
	 *
	 */
	template <typename id1, typename id2> inline const openfpm::vector<std::pair<size_t,size_t>> ghost_processorID_pair(Point<dim,T> & p, const int opt = MULTIPLE)
	{
		ids_p.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getCellIterator(geo_cell.getCell(p));

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInsideNP(p) == true)
			{
				ids_p.add(std::pair<size_t,size_t>(id1::id(vb_int.get(bid),bid),id2::id(vb_int.get(bid),bid)));
			}

			++cell_it;
		}

		// Make the id unique
		if (opt == UNIQUE)
		{
			ids_p.sort();
			ids_p.unique();
		}

		return ids_p;
	}

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 * (Internal ghost)
	 *
	 * if the particle come from an internal ghost from the periodicity of the domain, position must be shifted
	 * this function return the id of the shift vector
	 *
	 * \see getShiftVector
	 *
	 * \tparam id type of id to get box_id processor_id lc_processor_id shift_id
	 * \param p Particle position
	 * \param opt intersection boxes of the same processor can overlap, so in general the function
	 *        can produce more entry with the same processor, the UNIQUE option eliminate double entries
	 *        (UNIQUE) is for particle data (MULTIPLE) is for grid data [default MULTIPLE]
	 *
	 * \return the processor ids
	 *
	 */
	template <typename id> inline const openfpm::vector<size_t> ghost_processorID(const Point<dim,T> & p, const int opt = MULTIPLE)
	{
		ids.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getCellIterator(geo_cell.getCell(p));

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInsideNP(p) == true)
			{
				ids.add(id::id(vb_int.get(bid),bid));
			}

			++cell_it;
		}

		// Make the id unique
		if (opt == UNIQUE)
		{
			ids_p.sort();
			ids_p.unique();
		}

		return ids;
	}

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 * (Internal ghost)
	 *
	 * \tparam id1 first index type to get box_id processor_id lc_processor_id
	 * \tparam id2 second index type to get box_id processor_id lc_processor_id
	 *
	 * \param p Particle position
	 * \param opt indicate if the entries in the vector must be unique
	 *
	 * \return a vector of pair containing the requested information
	 *
	 */
	template<typename id1, typename id2, typename Mem> inline const openfpm::vector<std::pair<size_t,size_t>> & ghost_processorID_pair(const encapc<1,Point<dim,T>,Mem> & p, const int opt = MULTIPLE)
	{
		ids_p.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getCellIterator(geo_cell.getCell(p));

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInsideNP(p) == true)
			{
				ids_p.add(std::pair<size_t,size_t>(id1::id(vb_int.get(bid),bid),id2::id(vb_int.get(bid),bid)));
			}

			++cell_it;
		}

		// Make the id unique
		if (opt == UNIQUE)
		{
			ids_p.sort();
			ids_p.unique();
		}

		return ids_p;
	}

	/*! \brief Given a position it return if the position belong to any neighborhood processor ghost
	 * (Internal ghost)
	 *
	 * \tparam id type of if to get box_id processor_id lc_processor_id
	 *
	 * \param p Particle position
	 * \param opt it indicate if the entry in the vector must be unique or not
	 *
	 * \return the processor ids
	 *
	 */
	template<typename id, typename Mem> inline const openfpm::vector<size_t> & ghost_processorID(const encapc<1,Point<dim,T>,Mem> & p, const int opt = MULTIPLE)
	{
		ids.clear();

		// Check with geo-cell if a particle is inside one Cell containing boxes

		auto cell_it = geo_cell.getCellIterator(geo_cell.getCell(p));

		// For each element in the cell, check if the point is inside the box
		// if it is, store the processor id
		while (cell_it.isNext())
		{
			size_t bid = cell_it.get();

			if (vb_int.get(bid).box.isInsideNP(p) == true)
			{
				ids.add(id::id(vb_int.get(bid),bid));
			}

			++cell_it;
		}

		// Make the id unique
		if (opt == UNIQUE)
		{
			ids_p.sort();
			ids_p.unique();
		}

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
	 *
	 * \return true if the write succeed
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

	/*! \brief Check if the ie_ghosts contain the same information
	 *
	 * \param ig Element to check
	 *
	 * \return true if they are equal
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
			if (getIGhostBox(i) != ig.getIGhostBox(i))
				return false;
			if (getIGhostBoxProcessor(i) != ig.getIGhostBoxProcessor(i))
				return false;
		}

		for (size_t i = 0 ; i < proc_int_box.size() ; i++)
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
		}

		for (size_t i = 0 ; i < getNEGhostBox() ; i++)
		{
			if (getEGhostBox(i) != ig.getEGhostBox(i))
				return false;
			if (getEGhostBoxProcessor(i) != ig.getEGhostBoxProcessor(i))
				return false;
		}

		for (size_t i = 0 ; i < proc_int_box.size() ; i++)
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
		}

		return true;
	}

	/*! \brief Check if the ie_loc_ghosts contain the same information with the exception of the ghost part
	 * It is anyway required that the ghost come from the same sub-domains decomposition
	 *
	 * \param ig Element to check
	 *
	 * \return true if they are equal
	 *
	 */
	bool is_equal_ng(ie_ghost<dim,T> & ig)
	{
		return true;
	}

	/*! \brief Reset the nn_prcs structure
	 *
	 */
	void reset()
	{
		box_nn_processor_int.clear();
		proc_int_box.clear();
		vb_ext.clear();
		vb_int.clear();
		geo_cell.clear();
		shifts.clear();
		ids_p.clear();
		ids.clear();
	}
};


#endif /* SRC_DECOMPOSITION_IE_GHOST_HPP_ */
