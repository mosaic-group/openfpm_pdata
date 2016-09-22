/*
 * vector_dist_comm.hpp
 *
 *  Created on: Aug 18, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_COMM_HPP_
#define SRC_VECTOR_VECTOR_DIST_COMM_HPP_

#define V_SUB_UNIT_FACTOR 64

#define SKIP_LABELLING 512

#define NO_POSITION 1
#define WITH_POSITION 2

/*! \brief This class is an helper for the communication of vector_dist
 *
 * \tparam dim Dimensionality of the space where the elements lives
 * \tparam St type of space float, double ...
 * \tparam prop properties the vector element store in OpenFPM data structure format
 * \tparam Decomposition Decomposition strategy to use CartDecomposition ...
 * \tparam Memory Memory pool where store the information HeapMemory ...
 *
 * \see vector_dist
 *
 */

template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St>, typename Memory = HeapMemory>
class vector_dist_comm
{
	//! definition of the send vector for position
	typedef openfpm::vector<Point<dim, St>, Memory> send_pos_vector;

	//! VCluster
	Vcluster & v_cl;

	//! Domain decomposition
	Decomposition dec;

	//! It map the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! For each near processor, outgoing particle id
	//! \warning opart is assumed to be an ordered list
	//! first id particle id
	//! second id shift id
	//! third id is the processor id
	openfpm::vector<aggregate<size_t,size_t,size_t>> m_opart;

	//! Per processor ordered particles id for ghost_get (see prc_g_opart)
	//! For each processor the internal vector store the id of the
	//! particles that must be communicated to the other processors
	openfpm::vector<openfpm::vector<aggregate<size_t,size_t>>> g_opart;

	// processor rank list of g_opart
	openfpm::vector<size_t> prc_g_opart;

	//! Sending buffer for the ghost particles position
	openfpm::vector<send_pos_vector> g_pos_send;

	//! It store the list of processor that communicate with us (local processor)
	//! from the last ghost get
	openfpm::vector<size_t> prc_recv_get;

	//! the same as prc_recv_get but for put
	openfpm::vector<size_t> prc_recv_put;

	//! the same as prc_recv_get but for map
	openfpm::vector<size_t> prc_recv_map;

	//! It store the size of the elements added for each processor that communicate with us (local processor)
	//! from the last ghost get
	openfpm::vector<size_t> recv_sz_get;

	//! The same as recv_sz_get but for put
	openfpm::vector<size_t> recv_sz_put;

	//! The same as recv_sz_get but for map
	openfpm::vector<size_t> recv_sz_map;

	//! Local ghost marker (across the ghost particles it mark from where we have the)
	//! replicated ghost particles that are local
	size_t lg_m;

	//! process the particle with properties
	struct proc_without_prp
	{
		template<typename T1, typename T2> inline static void proc(size_t lbl, size_t cnt, size_t id, T1 & v_prp, T2 & m_prp)
		{
			m_prp.get(lbl).set(cnt, v_prp.get(id));
		}
	};

	template<typename prp_object, int ... prp>
	struct proc_with_prp
	{
		template<typename T1, typename T2> inline static void proc(size_t lbl, size_t cnt, size_t id, T1 & v_prp, T2 & m_prp)
		{
			// source object type
			typedef encapc<1, prop, typename openfpm::vector<prop>::layout_type> encap_src;
			// destination object type
			typedef encapc<1, prp_object, typename openfpm::vector<prp_object>::layout_type> encap_dst;

			// Copy only the selected properties
			object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(id), m_prp.get(lbl).get(cnt));
		}
	};

	//! It process one particle
	template<typename proc_class, typename T1, typename T2, typename T3, typename T4> inline void process_map_particle(size_t i, long int & end, long int & id_end, T1 & m_pos, T2 & m_prp, T3 & v_pos, T4 & v_prp, openfpm::vector<size_t> & cnt)
	{
		long int prc_id = m_opart.template get<2>(i);
		size_t id = m_opart.template get<0>(i);

		if (prc_id >= 0)
		{
			size_t lbl = p_map_req.get(prc_id);

			m_pos.get(lbl).set(cnt.get(lbl), v_pos.get(id));
			proc_class::proc(lbl,cnt.get(lbl),id,v_prp,m_prp);

			cnt.get(lbl)++;

			// swap the particle
			long int id_valid = get_end_valid(end,id_end);

			if (id_valid > 0 && (long int)id < id_valid)
			{
				v_pos.set(id,v_pos.get(id_valid));
				v_prp.set(id,v_prp.get(id_valid));
			}
		}
		else
		{
			// swap the particle
			long int id_valid = get_end_valid(end,id_end);

			if (id_valid > 0 && (long int)id < id_valid)
			{
				v_pos.set(id,v_pos.get(id_valid));
				v_prp.set(id,v_prp.get(id_valid));
			}
		}
	}

	/*! \brief Return a valid particle starting from end and tracing back
	 *
	 * \param end actual opart particle pointer
	 * \param actual end particle point
	 *
	 * \return a valid particle
	 *
	 */
	inline size_t get_end_valid(long int & end, long int & end_id)
	{
		end_id--;

		while (end >= 0 && end_id >= 0 && (long int)m_opart.template get<0>(end) == end_id)
		{
			end_id--;
			end--;
		}

		return end_id;
	}

	/*! \brief It store for each processor the position and properties vector of the particles
	 *
	 * This structure is used in the map function
	 *
	 */
	struct pos_prop
	{
		//! position vector
		openfpm::vector<Point<dim, St>, PreAllocHeapMemory<2>, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity> pos;
		//! properties vector
		openfpm::vector<prop, PreAllocHeapMemory<2>, typename memory_traits_lin<prop>::type, memory_traits_lin, openfpm::grow_policy_identity> prp;
	};

	/*! \brief for each processor store 2 vector containing the sending buffers
	 *
	 * This structure is used in the map_list function
	 *
	 */
	template <typename sel_prop>
	struct pos_prop_sel
	{
		//! position vector
		openfpm::vector<Point<dim, St>, PreAllocHeapMemory<2>, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity> pos;
		//! properties vector
		openfpm::vector<sel_prop, PreAllocHeapMemory<2>, typename memory_traits_lin<sel_prop>::type, memory_traits_lin, openfpm::grow_policy_identity> prp;
	};

	//! Flags that indicate that the function createShiftBox() has been called
	bool is_shift_box_created = false;

	//! this map is used to check if a combination is already present
	std::unordered_map<size_t, size_t> map_cmb;

	//! The boxes touching the border of the domain are divided in groups (first vector)
	//! each group contain internal ghost coming from sub-domains of the same section
	openfpm::vector_std<openfpm::vector_std<Box<dim, St>>>box_f;

	//! Store the sector for each group (previous vector)
	openfpm::vector_std<comb<dim>> box_cmb;

	//! Id of the local particle to replicate for ghost_get
	openfpm::vector<aggregate<size_t,size_t>> o_part_loc;

	/*! \brief For every internal ghost box we create a structure that order such internal local ghost box in
	 *         shift vectors
	 *
	 */
	void createShiftBox()
	{
		if (is_shift_box_created == true)
			return;

		// Add local particles coming from periodic boundary, the only boxes that count are the one
		// touching the border, filter them
		for (size_t i = 0; i < dec.getNLocalSub(); i++)
		{
			size_t Nl = dec.getLocalNIGhost(i);

			for (size_t j = 0; j < Nl; j++)
			{
				// If the ghost does not come from the intersection with an out of
				// border sub-domain the combination is all zero and n_zero return dim
				if (dec.getLocalIGhostPos(i, j).n_zero() == dim)
					continue;

				// Check if we already have boxes with such combination
				auto it = map_cmb.find(dec.getLocalIGhostPos(i, j).lin());
				if (it == map_cmb.end())
				{
					// we do not have it
					box_f.add();
					box_f.last().add(dec.getLocalIGhostBox(i, j));
					box_cmb.add(dec.getLocalIGhostPos(i, j));
					map_cmb[dec.getLocalIGhostPos(i, j).lin()] = box_f.size() - 1;
				}
				else
				{
					// we have it
					box_f.get(it->second).add(dec.getLocalIGhostBox(i, j));
				}

			}
		}

		is_shift_box_created = true;
	}

	/*! \brief Local ghost from labeled particles
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particles properties
	 *
	 */
	void local_ghost_from_opart(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp)
	{
		// get the shift vectors
		const openfpm::vector<Point<dim, St>> & shifts = dec.getShiftVectors();

		for (size_t i = 0 ; i < o_part_loc.size() ; i++)
		{
			size_t lin_id = o_part_loc.get<1>(i);
			size_t key = o_part_loc.template get<0>(i);

			Point<dim, St> p = v_pos.get(key);
			// shift
			p -= shifts.get(lin_id);

			// add this particle shifting its position
			v_pos.add(p);
			v_prp.add();
			v_prp.last() = v_prp.get(key);
		}
	}

	/*! \brief Local ghost from decomposition
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	void local_ghost_from_dec(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t g_m)
	{
		o_part_loc.clear();

		// get the shift vectors
		const openfpm::vector<Point<dim, St>> & shifts = dec.getShiftVectors();

		// Label the internal (assigned) particles
		auto it = v_pos.getIteratorTo(g_m);

		while (it.isNext())
		{
			auto key = it.get();

			// If particles are inside these boxes
			for (size_t i = 0; i < box_f.size(); i++)
			{
				for (size_t j = 0; j < box_f.get(i).size(); j++)
				{
					if (box_f.get(i).get(j).isInside(v_pos.get(key)) == true)
					{
						size_t lin_id = box_cmb.get(i).lin();

						o_part_loc.add();
						o_part_loc.template get<0>(o_part_loc.size()-1) = key;
						o_part_loc.template get<1>(o_part_loc.size()-1) = lin_id;

						Point<dim, St> p = v_pos.get(key);
						// shift
						p -= shifts.get(lin_id);

						// add this particle shifting its position
						v_pos.add(p);
						v_prp.add();
						v_prp.last() = v_prp.get(key);

						// boxes in one group can be overlapping
						// we do not have to search for the other
						// boxes otherwise we will have duplicate particles
						//
						// A small note overlap of boxes across groups is fine
						// (and needed) because each group has different shift
						// producing non overlapping particles
						//
						break;
					}
				}
			}

			++it;
		}
	}

	/*! \brief Add local particles based on the boundary conditions
	 *
	 * In order to understand what this function use the following
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
 B	 |         |                        |     A   |
 *	 |         |                        |    *    |
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
	 *  The box is the domain, while all boxes at the border (so not (0,0) ) are the
	 *  ghost part at the border of the domain. If a particle A is in the position in figure
	 *  a particle B must be created. This function duplicate the particle A, if A and B are
	 *  local
	 *
	 * \param v_pos vector of particle of positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 * \param opt options
	 *
	 */
	void add_loc_particles_bc(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp ,size_t & g_m, size_t opt)
	{
		// Create the shift boxes
		createShiftBox();

		lg_m = v_prp.size();

		if (box_f.size() == 0)
			return;
		else
		{
			if (opt & SKIP_LABELLING)
				local_ghost_from_opart(v_pos,v_prp);
			else
				local_ghost_from_dec(v_pos,v_prp,g_m);
		}
	}

	/*! \brief This function fill the send buffer for the particle position after the particles has been label with labelParticles
	 *
	 * \param v_pos vector of particle positions
	 * \param g_pos_send Send buffer to fill
	 * \param prAlloc_pos Memory object for the send buffer
	 *
	 */
	void fill_send_ghost_pos_buf(openfpm::vector<Point<dim, St>> & v_pos,openfpm::vector<send_pos_vector> & g_pos_send)
	{
		// get the shift vectors
		const openfpm::vector<Point<dim, St>> & shifts = dec.getShiftVectors();

		// create a number of send buffers equal to the near processors
		g_pos_send.resize(g_opart.size());
		for (size_t i = 0; i < g_pos_send.size(); i++)
		{
			// resize the sending vector (No allocation is produced)
			g_pos_send.get(i).resize(g_opart.get(i).size());
		}

		// Fill the send buffer
		for (size_t i = 0; i < g_opart.size(); i++)
		{
			for (size_t j = 0; j < g_opart.get(i).size(); j++)
			{
				Point<dim, St> s = v_pos.get(g_opart.get(i).template get<0>(j));
				s -= shifts.get(g_opart.get(i).template get<1>(j));
				g_pos_send.get(i).set(j, s);
			}
		}
	}

	/*! \brief This function fill the send buffer for ghost_put
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param v_prp vector of particle properties
	 * \param g_send_prp Send buffer to fill
	 * \param prAlloc_prp Memory object for the send buffer
	 * \param g_m ghost marker
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp> void fill_send_ghost_put_prp_buf(openfpm::vector<prop> & v_prp, openfpm::vector<send_vector> & g_send_prp, size_t & g_m)
	{
		// create a number of send buffers equal to the near processors
		// from which we received
		g_send_prp.resize(prc_recv_get.size());
		for (size_t i = 0; i < g_send_prp.size(); i++)
		{
			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(recv_sz_get.get(i));
		}

		size_t accum = g_m;

		// Fill the send buffer
		for (size_t i = 0; i < prc_recv_get.size(); i++)
		{
			size_t j2 = 0;
			for (size_t j = accum; j < accum + recv_sz_get.get(i); j++)
			{
				// source object type
				typedef encapc<1, prop, typename openfpm::vector<prop>::layout_type> encap_src;
				// destination object type
				typedef encapc<1, prp_object, typename openfpm::vector<prp_object>::layout_type> encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(j), g_send_prp.get(i).get(j2));

				j2++;
			}

			accum = accum + recv_sz_get.get(i);
		}
	}

	/*! \brief This function fill the send buffer for properties after the particles has been label with labelParticles
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param v_prp vector of particle properties
	 * \param g_send_prp Send buffer to fill
	 * \param prAlloc_prp Memory object for the send buffer
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp> void fill_send_ghost_prp_buf(openfpm::vector<prop> & v_prp, openfpm::vector<send_vector> & g_send_prp)
	{
		// create a number of send buffers equal to the near processors
		g_send_prp.resize(g_opart.size());
		for (size_t i = 0; i < g_send_prp.size(); i++)
		{
			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(g_opart.get(i).size());
		}

		// Fill the send buffer
		for (size_t i = 0; i < g_opart.size(); i++)
		{
			for (size_t j = 0; j < g_opart.get(i).size(); j++)
			{
				// source object type
				typedef encapc<1, prop, typename openfpm::vector<prop>::layout_type> encap_src;
				// destination object type
				typedef encapc<1, prp_object, typename openfpm::vector<prp_object>::layout_type> encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(g_opart.get(i).template get<0>(j)), g_send_prp.get(i).get(j));
			}
		}
	}

	/*! \brief allocate and fill the send buffer for the map function
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particles properties
	 * \param prc_r List of processor rank involved in the send
	 * \param prc_sz_r For each processor in the list the size of the message to send
	 * \param pb send buffer
	 *
	 */
	void fill_send_map_buf(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<openfpm::vector<Point<dim,St>>> & m_pos, openfpm::vector<openfpm::vector<prop>> & m_prp)
	{
		m_prp.resize(prc_sz_r.size());
		m_pos.resize(prc_sz_r.size());
		openfpm::vector<size_t> cnt(prc_sz_r.size());

		for (size_t i = 0; i < prc_sz_r.size() ; i++)
		{
			// set the size and allocate, using mem warant that pos and prp is contiguous
			m_pos.get(i).resize(prc_sz_r.get(i));
			m_prp.get(i).resize(prc_sz_r.get(i));
			cnt.get(i) = 0;
		}

		// end vector point
		long int id_end = v_pos.size();

		// end opart point
		long int end = m_opart.size()-1;

		// Run through all the particles and fill the sending buffer
		for (size_t i = 0; i < m_opart.size(); i++)
		{
			process_map_particle<proc_without_prp>(i,end,id_end,m_pos,m_prp,v_pos,v_prp,cnt);
		}

		v_pos.resize(v_pos.size() - m_opart.size());
		v_prp.resize(v_prp.size() - m_opart.size());
	}


	/*! \brief allocate and fill the send buffer for the map function
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param prc_r List of processor rank involved in the send
	 * \param prc_sz_r For each processor in the list the size of the message to send
	 * \param pb send buffer
	 *
	 */
	template<typename prp_object,int ... prp> void fill_send_map_buf_list(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<openfpm::vector<Point<dim,St>>> & m_pos, openfpm::vector<openfpm::vector<prp_object>> & m_prp)
	{
		m_prp.resize(prc_sz_r.size());
		m_pos.resize(prc_sz_r.size());
		openfpm::vector<size_t> cnt(prc_sz_r.size());

		for (size_t i = 0; i < prc_sz_r.size(); i++)
		{
			// set the size and allocate, using mem warant that pos and prp is contiguous
			m_pos.get(i).resize(prc_sz_r.get(i));
			m_prp.get(i).resize(prc_sz_r.get(i));
			cnt.get(i) = 0;
		}

		// end vector point
		long int id_end = v_pos.size();

		// end opart point
		long int end = m_opart.size()-1;

		// Run through all the particles and fill the sending buffer
		for (size_t i = 0; i < m_opart.size(); i++)
		{
			process_map_particle<proc_with_prp<prp_object,prp...>>(i,end,id_end,m_pos,m_prp,v_pos,v_prp,cnt);
		}

		v_pos.resize(v_pos.size() - m_opart.size());
		v_prp.resize(v_prp.size() - m_opart.size());
	}

	/*! \brief Label particles for mappings
	 *
	 * \param v_pos vector of particle positions
	 * \param lbl_p Particle labeled
	 * \param prc_sz For each processor the number of particles to send
	 * \param opart id of the particles to send
	 *
	 */
	template<typename obp> void labelParticleProcessor(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<aggregate<size_t,size_t,size_t>> & lbl_p, openfpm::vector<size_t> & prc_sz)
	{
		// reset lbl_p
		lbl_p.clear();

		// resize the label buffer
		prc_sz.resize(v_cl.getProcessingUnits());

		auto it = v_pos.getIterator();

		// Label all the particles with the processor id where they should go
		while (it.isNext())
		{
			auto key = it.get();

			// Apply the boundary conditions
			dec.applyPointBC(v_pos.get(key));

			size_t p_id = 0;

			// Check if the particle is inside the domain
			if (dec.getDomain().isInside(v_pos.get(key)) == true)
				p_id = dec.processorIDBC(v_pos.get(key));
			else
				p_id = obp::out(key, v_cl.getProcessUnitID());

			// Particle to move
			if (p_id != v_cl.getProcessUnitID())
			{
				if ((long int) p_id != -1)
				{
					prc_sz.get(p_id)++;
					lbl_p.add();
					lbl_p.last().template get<0>() = key;
					lbl_p.last().template get<2>() = p_id;
				}
				else
				{
					lbl_p.add();
					lbl_p.last().template get<0>() = key;
					lbl_p.last().template get<2>() = p_id;
				}
			}

			// Add processors and add size

			++it;
		}
	}

	/*! \brief Label the particles
	 *
	 * It count the number of particle to send to each processors and save its ids
	 *
	 * \see nn_prcs::getShiftvectors()
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	void labelParticlesGhost(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, openfpm::vector<size_t> & prc, size_t & g_m)
	{
		// Buffer that contain for each processor the id of the particle to send
		g_opart.clear();
		g_opart.resize(dec.getNNProcessors());
		prc_g_opart.clear();

		// Iterate over all particles
		auto it = v_pos.getIteratorTo(g_m);
		while (it.isNext())
		{
			auto key = it.get();

			// Given a particle, it return which processor require it (first id) and shift id, second id
			// For an explanation about shifts vectors please consult getShiftVector in ie_ghost
			const openfpm::vector<std::pair<size_t, size_t>> & vp_id = dec.template ghost_processorID_pair<typename Decomposition::lc_processor_id, typename Decomposition::shift_id>(v_pos.get(key), UNIQUE);

			for (size_t i = 0; i < vp_id.size(); i++)
			{
				// processor id
				size_t p_id = vp_id.get(i).first;

				// add particle to communicate
				g_opart.get(p_id).add();
				g_opart.get(p_id).last().template get<0>() = key;
				g_opart.get(p_id).last().template get<1>() = vp_id.get(i).second;
			}

			++it;
		}

		// remove all zero entry and construct prc (the list of the sending processors)
		openfpm::vector<openfpm::vector<aggregate<size_t,size_t>>> g_opart_f;

		// count the non zero element
		for (size_t i = 0 ; i < g_opart.size() ; i++)
		{
			if (g_opart.get(i).size() != 0)
			{
				g_opart_f.add();
				g_opart.get(i).swap(g_opart_f.last());
				prc.add(dec.IDtoProc(i));
			}
		}

		g_opart.swap(g_opart_f);
	}

	/*! \brief Call-back to allocate buffer to receive incoming elements (particles)
	 *
	 * \param msg_i size required to receive the message from i
	 * \param total_msg total size to receive from all the processors
	 * \param total_p the total number of processor that want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr a pointer to the vector_dist structure
	 *
	 * \return the pointer where to store the message for the processor i
	 *
	 */
	static void * message_alloc_map(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		// cast the pointer
		vector_dist_comm<dim, St, prop, Decomposition, Memory> * vd = static_cast<vector_dist_comm<dim, St, prop, Decomposition, Memory> *>(ptr);

		vd->recv_mem_gm.resize(vd->v_cl.getProcessingUnits());
		vd->recv_mem_gm.get(i).resize(msg_i);

		return vd->recv_mem_gm.get(i).getPointer();
	}

public:

	/*! \brief Copy Constructor
	 *
	 * \param v vector to copy
	 *
	 */
	vector_dist_comm(const vector_dist_comm<dim,St,prop,Decomposition,Memory> & v)
	:v_cl(create_vcluster()),dec(create_vcluster())
	{
		this->operator=(v);
	}


	/*! \brief Constructor
	 *
	 * \param dec Domain decompositon
	 *
	 */
	vector_dist_comm(const Decomposition & dec)
	:v_cl(create_vcluster()),dec(dec)
	{

	}

	/*! \brief Constructor
	 *
	 * \param dec Domain decompositon
	 *
	 */
	vector_dist_comm(Decomposition && dec)
	:v_cl(create_vcluster()),dec(dec)
	{

	}

	/*! \brief Constructor
	 *
	 */
	vector_dist_comm()
	:v_cl(create_vcluster()),dec(create_vcluster())
	{
	}

	/*! \brief Get the number of minimum sub-domain
	 *
	 * \return minimum number
	 *
	 */
	static size_t getDefaultNsubsub()
	{
		return V_SUB_UNIT_FACTOR;
	}

	/*! \brief Initialize the decomposition
	 *
	 * \param box domain
	 * \param bc boundary conditions
	 * \param g ghost extension
	 *
	 */
	void init_decomposition(Box<dim,St> & box, const size_t (& bc)[dim],const Ghost<dim,St> & g)
	{
		// Create a valid decomposition of the space
		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * getDefaultNsubsub();

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[dim];
		for (size_t i = 0; i < dim; i++)
		{
			div[i] = openfpm::math::round_big_2(pow(n_sub, 1.0 / dim));
		}

		// Create the sub-domains
		dec.setParameters(div, box, bc, g);
		dec.decompose();
	}

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam prp list of properties to get synchronize
	 *
	 * \param opt options WITH_POSITION, it send also the positional information of the particles
	 * \param v_pos vector of position to update
	 * \param v_prp vector of properties to update
	 * \param g_m marker between real and ghost particles
	 *
	 */
	template<int ... prp> inline void ghost_get_(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & g_m, size_t opt = WITH_POSITION)
	{
		// Unload receive buffer
		for (size_t i = 0 ; i < recv_sz_get.size() ; i++)
			recv_sz_get.get(i) = 0;

		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object> send_vector;

		// reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Label all the particles
		if ((opt & SKIP_LABELLING) == false)
			labelParticlesGhost(v_pos,v_prp,prc_g_opart,g_m);

		// Send and receive ghost particle information
		openfpm::vector<send_vector> g_send_prp;
		fill_send_ghost_prp_buf<send_vector, prp_object, prp...>(v_prp,g_send_prp);

		// Create and fill the send buffer for the particle position
		if (opt != NO_POSITION)
			fill_send_ghost_pos_buf(v_pos,g_pos_send);

		prc_recv_get.clear();
		recv_sz_get.clear();
		v_cl.SSendRecvP<send_vector,decltype(v_prp),prp...>(g_send_prp,v_prp,prc_g_opart,prc_recv_get,recv_sz_get);

		if (opt != NO_POSITION)
		{
			prc_recv_get.clear();
			recv_sz_get.clear();
			v_cl.SSendRecv(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get);
		}

		add_loc_particles_bc(v_pos,v_prp,g_m,opt);
	}


	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 * \tparam prp properties to communicate
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	template<unsigned int ... prp> void map_list_(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & g_m)
	{
		typedef KillParticle obp;

		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// It contain the list of the processors this processor should to communicate with
		openfpm::vector<size_t> p_list;

		// map completely reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Contain the processor id of each particle (basically where they have to go)
		labelParticleProcessor<obp>(v_pos,m_opart, prc_sz);

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		p_map_req.resize(v_cl.getProcessingUnits());
		openfpm::vector<size_t> prc_sz_r;
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (prc_sz.get(i) != 0)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(prc_sz.get(i));
			}
		}

		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		//! position vector
		openfpm::vector<openfpm::vector<Point<dim, St>>> m_pos;
		//! properties vector
		openfpm::vector<openfpm::vector<prp_object>> m_prp;

		fill_send_map_buf_list<prp_object,prp...>(v_pos,v_prp,prc_sz_r, m_pos, m_prp);

		v_cl.SSendRecv(m_pos,v_pos,prc_r,prc_recv_map,recv_sz_map);
		v_cl.SSendRecvP<openfpm::vector<prp_object>,decltype(v_prp),prp...>(m_prp,v_prp,prc_r,prc_recv_map,recv_sz_map);

		// mark the ghost part

		g_m = v_pos.size();
	}

	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	template<typename obp = KillParticle> void map_(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & g_m)
	{
		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// It contain the list of the processors this processor should to communicate with
		openfpm::vector<size_t> p_list;

		// map completely reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Contain the processor id of each particle (basically where they have to go)
		labelParticleProcessor<obp>(v_pos,m_opart, prc_sz);

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		p_map_req.resize(v_cl.getProcessingUnits());
		openfpm::vector<size_t> prc_sz_r;
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (prc_sz.get(i) != 0)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(prc_sz.get(i));
			}
		}

		//! position vector
		openfpm::vector<openfpm::vector<Point<dim, St>>> m_pos;
		//! properties vector
		openfpm::vector<openfpm::vector<prop>> m_prp;

		fill_send_map_buf(v_pos,v_prp, prc_sz_r, m_pos, m_prp);

		v_cl.SSendRecv(m_pos,v_pos,prc_r,prc_recv_map,recv_sz_map);
		v_cl.SSendRecv(m_prp,v_prp,prc_r,prc_recv_map,recv_sz_map);

		// mark the ghost part

		g_m = v_pos.size();
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	inline Decomposition & getDecomposition()
	{
		return dec;
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	inline const Decomposition & getDecomposition() const
	{
		return dec;
	}

	/*! \brief Copy a vector
	 *
	 * \param vc vector to copy
	 *
	 * \return iteself
	 *
	 */
	vector_dist_comm<dim,St,prop,Decomposition,Memory> & operator=(const vector_dist_comm<dim,St,prop,Decomposition,Memory> & vc)
	{
		dec = vc.dec;

		return *this;
	}

	/*! \brief Copy a vector
	 *
	 * \param vc vector to copy
	 *
	 * \return itself
	 *
	 */
	vector_dist_comm<dim,St,prop,Decomposition,Memory> & operator=(vector_dist_comm<dim,St,prop,Decomposition,Memory> && vc)
	{
		dec = vc.dec;

		return *this;
	}

	/*! \brief Ghost put
	 *
	 * \tparam op operation to apply
	 * \tparam prp set of properties
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector od particle properties
	 * \param g_m ghost marker
	 *
	 */
	template<template<typename,typename> class op, int ... prp> void ghost_put_(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & g_m)
	{
		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object> send_vector;

		openfpm::vector<send_vector> g_send_prp;
		fill_send_ghost_put_prp_buf<send_vector, prp_object, prp...>(v_prp,g_send_prp,g_m);

		// Send and receive ghost particle information
		op_ssend_recv_merge<op> opm(g_opart);
		v_cl.SSendRecvP_op<op_ssend_recv_merge<op>,send_vector,decltype(v_prp),prp...>(g_send_prp,v_prp,prc_recv_get,opm,prc_recv_put,recv_sz_put);
	}
};


#endif /* SRC_VECTOR_VECTOR_DIST_COMM_HPP_ */
