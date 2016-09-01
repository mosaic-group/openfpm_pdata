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
	//! VCluster
	Vcluster & v_cl;

	//! Domain decomposition
	Decomposition dec;

	//! It map the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! For each near processor, outgoing particle id
	openfpm::vector<openfpm::vector<size_t>> opart;

	//! For each near processor, particle shift vector
	openfpm::vector<openfpm::vector<size_t>> oshift;

	//! For each adjacent processor the size of the ghost sending buffer
	openfpm::vector<size_t> ghost_prc_sz;

	//! Sending buffer for the ghost particles properties
	BHeapMemory g_prp_mem;

	//! Sending buffer for the ghost particles position
	BHeapMemory g_pos_mem;

	//! For each adjacent processor it store from which processor come from
	openfpm::vector<size_t> prc_recv;

	//! the same as prc_recv but for put
	openfpm::vector<size_t> prc_recv_put;

	//! Number of received elements
	openfpm::vector<size_t> n_recv_ele;

	//! For each adjacent processor it store the size of the receiving message in byte
	openfpm::vector<size_t> recv_sz;

	//! The same as recv_sz but for put
	openfpm::vector<size_t> recv_sz_put;

	//! For each adjacent processor it store the received message for ghost get
	openfpm::vector<BHeapMemory> recv_mem_gg;

	//! For each processor it store the received message for global map
	openfpm::vector<BHeapMemory> recv_mem_gm;

	//! Local ghost marker (across the ghost particles it mark from where we have the)
	//! replicated ghost particles that are local
	size_t lg_m;

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

	//! definition of the send vector for position
	typedef openfpm::vector<Point<dim, St>, ExtPreAlloc<Memory>, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin , openfpm::grow_policy_identity> send_pos_vector;

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
	void fill_send_ghost_pos_buf(openfpm::vector<Point<dim, St>> & v_pos,openfpm::vector<send_pos_vector> & g_pos_send, ExtPreAlloc<Memory> * prAlloc_pos)
	{
		// get the shift vectors
		const openfpm::vector<Point<dim, St>> & shifts = dec.getShiftVectors();

		// create a number of send buffers equal to the near processors
		g_pos_send.resize(ghost_prc_sz.size());
		for (size_t i = 0; i < g_pos_send.size(); i++)
		{
			// set the preallocated memory to ensure contiguity
			g_pos_send.get(i).setMemory(*prAlloc_pos);

			// resize the sending vector (No allocation is produced)
			g_pos_send.get(i).resize(ghost_prc_sz.get(i));
		}

		// Fill the send buffer
		for (size_t i = 0; i < opart.size(); i++)
		{
			for (size_t j = 0; j < opart.get(i).size(); j++)
			{
				Point<dim, St> s = v_pos.get(opart.get(i).get(j));
				s -= shifts.get(oshift.get(i).get(j));
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
	template<typename send_vector, typename prp_object, int ... prp> void fill_send_ghost_put_prp_buf(openfpm::vector<prop> & v_prp, openfpm::vector<send_vector> & g_send_prp, ExtPreAlloc<Memory> * prAlloc_prp, size_t & g_m)
	{
		// create a number of send buffers equal to the near processors
		// from which we received
		g_send_prp.resize(prc_recv.size());
		for (size_t i = 0; i < g_send_prp.size(); i++)
		{
			// set the preallocated memory to ensure contiguity
			g_send_prp.get(i).setMemory(*prAlloc_prp);

			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(n_recv_ele.get(i));
		}

		size_t accum = g_m;

		// Fill the send buffer
		for (size_t i = 0; i < prc_recv.size(); i++)
		{
			size_t j2 = 0;
			for (size_t j = accum; j < accum + n_recv_ele.get(i); j++)
			{
				// source object type
				typedef encapc<1, prop, typename openfpm::vector<prop>::layout_type> encap_src;
				// destination object type
				typedef encapc<1, prp_object, typename openfpm::vector<prp_object>::layout_type> encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(j), g_send_prp.get(i).get(j2));

				j2++;
			}

			accum = accum + n_recv_ele.get(i);
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
	template<typename send_vector, typename prp_object, int ... prp> void fill_send_ghost_prp_buf(openfpm::vector<prop> & v_prp, openfpm::vector<send_vector> & g_send_prp, ExtPreAlloc<Memory> * prAlloc_prp)
	{
		// create a number of send buffers equal to the near processors
		g_send_prp.resize(ghost_prc_sz.size());
		for (size_t i = 0; i < g_send_prp.size(); i++)
		{
			// set the preallocated memory to ensure contiguity
			g_send_prp.get(i).setMemory(*prAlloc_prp);

			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(ghost_prc_sz.get(i));
		}

		// Fill the send buffer
		for (size_t i = 0; i < opart.size(); i++)
		{
			for (size_t j = 0; j < opart.get(i).size(); j++)
			{
				// source object type
				typedef encapc<1, prop, typename openfpm::vector<prop>::layout_type> encap_src;
				// destination object type
				typedef encapc<1, prp_object, typename openfpm::vector<prp_object>::layout_type> encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(opart.get(i).get(j)), g_send_prp.get(i).get(j));
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
	void fill_send_map_buf(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp,openfpm::vector<size_t> & prc_r, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<pos_prop> & pb)
	{
		pb.resize(prc_r.size());

		for (size_t i = 0; i < prc_r.size(); i++)
		{
			// Create the size required to store the particles position and properties to communicate
			size_t s1 = openfpm::vector<Point<dim, St>, HeapMemory, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i), 0);
			size_t s2 = openfpm::vector<prop, HeapMemory, typename memory_traits_lin<prop>::type, memory_traits_lin, openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i), 0);

			// Preallocate the memory
			size_t sz[2] = { s1, s2 };
			PreAllocHeapMemory<2> * mem = new PreAllocHeapMemory<2>(sz);

			// Set the memory allocator
			pb.get(i).pos.setMemory(*mem);
			pb.get(i).prp.setMemory(*mem);

			// set the size and allocate, using mem warant that pos and prp is contiguous
			pb.get(i).pos.resize(prc_sz_r.get(i));
			pb.get(i).prp.resize(prc_sz_r.get(i));
		}

		// Run through all the particles and fill the sending buffer

		for (size_t i = 0; i < opart.size(); i++)
		{
			auto it = opart.get(i).getIterator();
			size_t lbl = p_map_req.get(i);

			while (it.isNext())
			{
				size_t key = it.get();
				size_t id = opart.get(i).get(key);

				pb.get(lbl).pos.set(key, v_pos.get(id));
				pb.get(lbl).prp.set(key, v_prp.get(id));

				++it;
			}
		}
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
	template<typename prp_object,int ... prp> void fill_send_map_buf_list(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, openfpm::vector<size_t> & prc_r, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<pos_prop_sel<prp_object>> & pb)
	{
		pb.resize(prc_r.size());

		for (size_t i = 0; i < prc_r.size(); i++)
		{
			// Create the size required to store the particles position and properties to communicate
			size_t s1 = openfpm::vector<Point<dim, St>, HeapMemory, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i), 0);
			size_t s2 = openfpm::vector<prp_object, HeapMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin, openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i), 0);

			// Preallocate the memory
			size_t sz[2] = { s1, s2 };
			PreAllocHeapMemory<2> * mem = new PreAllocHeapMemory<2>(sz);

			// Set the memory allocator
			pb.get(i).pos.setMemory(*mem);
			pb.get(i).prp.setMemory(*mem);

			// set the size and allocate, using mem warant that pos and prp is contiguous
			pb.get(i).pos.resize(prc_sz_r.get(i));
			pb.get(i).prp.resize(prc_sz_r.get(i));
		}

		// Run through all the particles and fill the sending buffer

		for (size_t i = 0; i < opart.size(); i++)
		{
			auto it = opart.get(i).getIterator();
			size_t lbl = p_map_req.get(i);

			while (it.isNext())
			{
				size_t key = it.get();
				size_t id = opart.get(i).get(key);

				pb.get(lbl).pos.set(key, v_pos.get(id));

				// source object type
				typedef encapc<1, prop, typename openfpm::vector<prop>::layout_type> encap_src;
				// destination object type
				typedef encapc<1, prp_object, typename openfpm::vector<prp_object>::layout_type> encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(id), pb.get(lbl).prp.get(key));

				++it;
			}
		}
	}

	/*! \brief Label particles for mappings
	 *
	 * \param v_pos vector of particle positions
	 * \param lbl_p Particle labeled
	 * \param prc_sz For each processor the number of particles to send
	 * \param opart id of the particles to send
	 *
	 */
	template<typename obp> void labelParticleProcessor(openfpm::vector<Point<dim, St>> & v_pos,openfpm::vector<openfpm::vector<size_t>> & lbl_p, openfpm::vector<size_t> & prc_sz, openfpm::vector<size_t> & opart)
	{
		// reset lbl_p
		lbl_p.resize(v_cl.getProcessingUnits());
		for (size_t i = 0; i < lbl_p.size(); i++)
			lbl_p.get(i).clear();

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
					lbl_p.get(p_id).add(key);
					opart.add(key);
				}
				else
				{
					opart.add(key);
				}
			}

			// Add processors and add size

			++it;
		}
	}

	/*! \brief This function process the received data for the properties and populate the ghost
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp> void process_received_ghost_prp(openfpm::vector<prop> & v_prp, size_t & g_m)
	{
		n_recv_ele.resize(recv_mem_gg.size());

		// Mark the ghost part
		g_m = v_prp.size();

		// Process the received data (recv_mem_gg != 0 if you have data)
		for (size_t i = 0; i < dec.getNNProcessors() && recv_mem_gg.size() != 0; i++)
		{
			// calculate the number of received elements
			size_t n_ele = recv_sz.get(i) / sizeof(prp_object);

			// add the received particles to the vector
			PtrMemory * ptr1 = new PtrMemory(recv_mem_gg.get(i).getPointer(), recv_sz.get(i));

			// create vector representation to a piece of memory already allocated
			openfpm::vector<prp_object, PtrMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin , openfpm::grow_policy_identity> v2;

			v2.setMemory(*ptr1);

			// resize with the number of elements and store the number
			// or received elements
			v2.resize(n_ele);
			n_recv_ele.get(i) = n_ele;

			// Add the ghost particle
			v_prp.template add_prp<prp_object, PtrMemory, openfpm::grow_policy_identity, prp...>(v2);
		}
	}


	/*! \brief This function process the received data from ghost put
	 *
	 * \tparam op operation to do
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	template<template<typename,typename> class op, typename send_vector, typename prp_object, int ... prp> void process_received_put_ghost_prp(openfpm::vector<prop> & v_prp, size_t g_m)
	{
		// Process the received data (recv_mem_gg != 0 if you have data)
		for (size_t i = 0; i <  recv_sz_put.size(); i++)
		{
			// calculate the number of received elements
			size_t n_ele = recv_sz_put.get(i) / sizeof(prp_object);

			// add the received particles to the vector
			PtrMemory * ptr1 = new PtrMemory(recv_mem_gg.get(i).getPointer(), recv_sz_put.get(i));

			// create vector representation to a piece of memory already allocated
			openfpm::vector<prp_object, PtrMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin , openfpm::grow_policy_identity> v2;

			v2.setMemory(*ptr1);

			// resize with the number of elements
			v2.resize(n_ele);

			// Add the ghost particle
			v_prp.template merge_prp<op,prp_object, PtrMemory, openfpm::grow_policy_identity, prp...>(v2,opart.get(i));
		}

		// process also the local replicated particles

		size_t i2 = 0;

#ifdef SE_CLASS1

		if (v_prp.size() - lg_m != o_part_loc.size())
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " v_prp.size() - lg_m = " << v_prp.size() - lg_m << " != " << o_part_loc.size() << std::endl;

#endif

		for (size_t i = lg_m ; i < v_prp.size() ; i++)
		{
			auto dst = v_prp.get(o_part_loc.template get<0>(i2));
			auto src = v_prp.get(i);
			copy_cpu_encap_encap_op_prp<op,decltype(v_prp.get(0)),decltype(v_prp.get(0)),prp...> cp(src,dst);

			boost::mpl::for_each_ref< boost::mpl::range_c<int,0,sizeof...(prp)> >(cp);

			i2++;
		}
	}

	/*! \brief This function process the received data for the properties and populate the ghost
	 *
	 * \param v_pos vector of particle positions
	 *
	 */
	void process_received_ghost_pos(openfpm::vector<Point<dim, St>> & v_pos)
	{
		// Process the received data (recv_mem_gg != 0 if you have data)
		for (size_t i = 0; i < dec.getNNProcessors() && recv_mem_gg.size() != 0; i++)
		{
			// calculate the number of received elements
			size_t n_ele = recv_sz.get(i) / sizeof(Point<dim, St> );

			// add the received particles to the vector
			PtrMemory * ptr1 = new PtrMemory(recv_mem_gg.get(i).getPointer(), recv_sz.get(i));

			// create vector representation to a piece of memory already allocated

			openfpm::vector<Point<dim, St>, PtrMemory, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin , openfpm::grow_policy_identity> v2;

			v2.setMemory(*ptr1);

			// resize with the number of elements
			v2.resize(n_ele);

			// Add the ghost particle
			v_pos.template add<PtrMemory, openfpm::grow_policy_identity>(v2);
		}
	}

	/*! \brief Process the received particles
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param out_part list of the out-going particles
	 *
	 */
	void process_received_map(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp,openfpm::vector<size_t> & out_part)
	{
		size_t o_p_id = 0;

		for (size_t i = 0; i < recv_mem_gm.size(); i++)
		{
			// Get the number of elements

			size_t n_ele = recv_mem_gm.get(i).size() / (sizeof(Point<dim, St> ) + sizeof(prop));

			// Pointer of the received positions for each near processor
			void * ptr_pos = (unsigned char *) recv_mem_gm.get(i).getPointer();
			// Pointer of the received properties for each near processor
			void * ptr_prp = (unsigned char *) recv_mem_gm.get(i).getPointer() + n_ele * sizeof(Point<dim, St> );

			PtrMemory * ptr1 = new PtrMemory(ptr_pos, n_ele * sizeof(Point<dim, St> ));
			PtrMemory * ptr2 = new PtrMemory(ptr_prp, n_ele * sizeof(prop));

			// create vector representation to a piece of memory already allocated

			openfpm::vector<Point<dim, St>, PtrMemory, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin ,openfpm::grow_policy_identity> vpos;
			openfpm::vector<prop, PtrMemory, typename memory_traits_lin<prop>::type, memory_traits_lin ,openfpm::grow_policy_identity> vprp;

			vpos.setMemory(*ptr1);
			vprp.setMemory(*ptr2);

			vpos.resize(n_ele);
			vprp.resize(n_ele);

			// Add the received particles to v_pos and v_prp

			size_t j = 0;
			for (; j < vpos.size() && o_p_id < out_part.size(); j++, o_p_id++)
			{
				v_pos.set(out_part.get(o_p_id), vpos.get(j));
				v_prp.set(out_part.get(o_p_id), vprp.get(j));
			}

			for (; j < vpos.size(); j++)
			{
				v_pos.add();
				v_pos.set(v_pos.size() - 1, vpos.get(j));
				v_prp.add();
				v_prp.set(v_prp.size() - 1, vprp.get(j));
			}
		}

		// remove the (out-going particles) in the vector

		v_pos.remove(out_part, o_p_id);
		v_prp.remove(out_part, o_p_id);
	}

	/*! \brief Process the received particles
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param out_part list of the out-going particles
	 *
	 */
	template<typename prp_object , int ... prp> void process_received_map_list(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, openfpm::vector<size_t> & out_part)
	{
		size_t o_p_id = 0;

		for (size_t i = 0; i < recv_mem_gm.size(); i++)
		{
			// Get the number of elements

			size_t n_ele = recv_mem_gm.get(i).size() / (sizeof(Point<dim, St> ) + sizeof(prp_object));

			// Pointer of the received positions for each near processor
			void * ptr_pos = (unsigned char *) recv_mem_gm.get(i).getPointer();
			// Pointer of the received properties for each near processor
			void * ptr_prp = (unsigned char *) recv_mem_gm.get(i).getPointer() + n_ele * sizeof(Point<dim, St> );

			PtrMemory * ptr1 = new PtrMemory(ptr_pos, n_ele * sizeof(Point<dim, St> ));
			PtrMemory * ptr2 = new PtrMemory(ptr_prp, n_ele * sizeof(prp_object));

			// create vector representation to a piece of memory already allocated

			openfpm::vector<Point<dim, St>, PtrMemory, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin ,openfpm::grow_policy_identity> vpos;
			openfpm::vector<prp_object, PtrMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin ,openfpm::grow_policy_identity> vprp;

			vpos.setMemory(*ptr1);
			vprp.setMemory(*ptr2);

			vpos.resize(n_ele);
			vprp.resize(n_ele);

			// Add the received particles to v_pos and v_prp

			size_t j = 0;
			for (; j < vpos.size() && o_p_id < out_part.size(); j++, o_p_id++)
			{
				v_pos.set(out_part.get(o_p_id), vpos.get(j));
				v_prp.template set_o<decltype(vprp.get(j)), prp... >(out_part.get(o_p_id), vprp.get(j));
			}

			for (; j < vpos.size(); j++)
			{
				v_pos.add();
				v_pos.set(v_pos.size() - 1, vpos.get(j));
				v_prp.template set_o<decltype(vprp.get(j)), prp... >(v_prp.size() - 1, vprp.get(j));
			}
		}

		// remove the (out-going particles) in the vector

		v_pos.remove(out_part, o_p_id);
		v_prp.remove(out_part, o_p_id);
	}

	/*! \brief Calculate send buffers total size and allocation
	 *
	 * \tparam prp_object object containing only the properties to send
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param size_byte_prp total size for the property buffer
	 * \param size_byte_pos total size for the position buffer
	 *
	 */
	template<typename prp_object> void calc_send_ghost_buf(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & size_byte_prp, size_t & size_byte_pos)
	{
		// Calculate the total size required for the sending buffer
		for (size_t i = 0; i < ghost_prc_sz.size(); i++)
		{
			size_t alloc_ele = openfpm::vector<prp_object, HeapMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin , openfpm::grow_policy_identity>::calculateMem(ghost_prc_sz.get(i), 0);
			size_byte_prp += alloc_ele;

			alloc_ele = openfpm::vector<Point<dim, St>, HeapMemory, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity>::calculateMem(ghost_prc_sz.get(i), 0);
			size_byte_pos += alloc_ele;
		}
	}

	/*! \brief Calculate send buffers total size and allocation
	 *
	 * \tparam prp_object object containing only the properties to send
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param size_byte_prp total size for the property buffer
	 * \param size_byte_pos total size for the position buffer
	 *
	 */
	template<typename prp_object> void calc_send_ghost_put_buf(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & size_byte_prp, size_t & size_byte_pos)
	{
		// Calculate the total size required for the sending buffer
		for (size_t i = 0; i < recv_sz.size(); i++)
		{
			size_t alloc_ele = openfpm::vector<prp_object, HeapMemory, typename memory_traits_lin<prp_object>::type, memory_traits_lin , openfpm::grow_policy_identity>::calculateMem(n_recv_ele.get(i), 0);
			size_byte_prp += alloc_ele;
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
	void labelParticlesGhost(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & g_m)
	{
		// Buffer that contain the number of elements to send for each processor
		ghost_prc_sz.clear();
		ghost_prc_sz.resize(dec.getNNProcessors());

		// Buffer that contain for each processor the id of the particle to send
		opart.clear();
		opart.resize(dec.getNNProcessors());

		// Buffer that contain for each processor the id of the shift vector
		oshift.clear();
		oshift.resize(dec.getNNProcessors());

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
				ghost_prc_sz.get(p_id)++;
				opart.get(p_id).add(key);
				oshift.get(p_id).add(vp_id.get(i).second);
			}

			++it;
		}
	}

	/*! \brief Call-back to allocate buffer to receive incoming elements (particles)
	 *
	 * \param msg_i message size required to receive from i
	 * \param total_msg message size to receive from all the processors
	 * \param total_p the total number of processor want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr void pointer parameter for additional data to pass to the call-back
	 *
	 */
	static void * msg_alloc_ghost_put(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		vector_dist_comm<dim, St, prop, Decomposition, Memory> * v = static_cast<vector_dist_comm<dim, St, prop, Decomposition, Memory> *>(ptr);

		v->recv_sz_put.resize(v->dec.getNNProcessors());
		v->recv_mem_gg.resize(v->dec.getNNProcessors());
		v->prc_recv_put.resize(v->dec.getNNProcessors());

		// Get the local processor id
		size_t lc_id = v->dec.ProctoID(i);

		// resize the receive buffer
		v->recv_mem_gg.get(lc_id).resize(msg_i);
		v->recv_sz_put.get(lc_id) = msg_i;

		// save the processor id
		v->prc_recv_put.get(lc_id) = i;

		return v->recv_mem_gg.get(lc_id).getPointer();
	}

	/*! \brief Call-back to allocate buffer to receive incoming elements (particles)
	 *
	 * \param msg_i message size required to receive from i
	 * \param total_msg message size to receive from all the processors
	 * \param total_p the total number of processor want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr void pointer parameter for additional data to pass to the call-back
	 *
	 */
	static void * msg_alloc_ghost_get(size_t msg_i, size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		vector_dist_comm<dim, St, prop, Decomposition, Memory> * v = static_cast<vector_dist_comm<dim, St, prop, Decomposition, Memory> *>(ptr);

		v->recv_sz.resize(v->dec.getNNProcessors());
		v->recv_mem_gg.resize(v->dec.getNNProcessors());
		v->prc_recv.resize(v->dec.getNNProcessors());

		// Get the local processor id
		size_t lc_id = v->dec.ProctoID(i);

		// resize the receive buffer
		v->recv_mem_gg.get(lc_id).resize(msg_i);
		v->recv_sz.get(lc_id) = msg_i;

		// save the processor id
		v->prc_recv.get(lc_id) = i;

		return v->recv_mem_gg.get(lc_id).getPointer();
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
		for (size_t i = 0 ; i < recv_sz.size() ; i++)
			recv_sz.get(i) = 0;

		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object, ExtPreAlloc<Memory>, typename memory_traits_lin<prp_object>::type, memory_traits_lin, openfpm::grow_policy_identity> send_vector;

		// reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Label all the particles
		if ((opt & SKIP_LABELLING) == false)
			labelParticlesGhost(v_pos,v_prp,g_m);

		// Calculate memory and allocation for the send buffers

		// Total size
		size_t size_byte_prp = 0;
		size_t size_byte_pos = 0;

		calc_send_ghost_buf<prp_object>(v_pos,v_prp,size_byte_prp, size_byte_pos);

		// Create memory for the send buffer

		g_prp_mem.resize(size_byte_prp);
		if (opt != NO_POSITION)
			g_pos_mem.resize(size_byte_pos);

		// Create and fill send buffer for particle properties

		ExtPreAlloc<Memory> * prAlloc_prp = new ExtPreAlloc<Memory>(size_byte_prp, g_prp_mem);

		openfpm::vector<send_vector> g_send_prp;
		fill_send_ghost_prp_buf<send_vector, prp_object, prp...>(v_prp,g_send_prp, prAlloc_prp);

		// Create and fill the send buffer for the particle position

		ExtPreAlloc<Memory> * prAlloc_pos;
		openfpm::vector<send_pos_vector> g_pos_send;
		if (opt != NO_POSITION)
		{
			prAlloc_pos = new ExtPreAlloc<Memory>(size_byte_pos, g_pos_mem);
			fill_send_ghost_pos_buf(v_pos,g_pos_send, prAlloc_pos);
		}

		// Create processor list
		openfpm::vector<size_t> prc;
		for (size_t i = 0; i < opart.size(); i++)
			prc.add(dec.IDtoProc(i));

		// Send/receive the particle properties information
		v_cl.sendrecvMultipleMessagesNBX(prc, g_send_prp, msg_alloc_ghost_get, this);
		process_received_ghost_prp<send_vector, prp_object, prp...>(v_prp,g_m);

		if (opt != NO_POSITION)
		{
			// Send/receive the particle properties information
			v_cl.sendrecvMultipleMessagesNBX(prc, g_pos_send, msg_alloc_ghost_get, this);
			process_received_ghost_pos(v_pos);
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

		// outgoing particles-id
		openfpm::vector<size_t> out_part;

		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// It contain the list of the processors this processor should to communicate with
		openfpm::vector<size_t> p_list;

		// map completely reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Contain the processor id of each particle (basically where they have to go)
		labelParticleProcessor<obp>(v_pos,opart, prc_sz, out_part);

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

		// Allocate the send buffers

		openfpm::vector<pos_prop_sel<prp_object>> pb;

		// fill the send buffers
		fill_send_map_buf_list<prp_object,prp...>(v_pos,v_prp,prc_r, prc_sz_r, pb);

		// Create the set of pointers
		openfpm::vector<void *> ptr(prc_r.size());
		for (size_t i = 0; i < prc_r.size(); i++)
		{
			ptr.get(i) = pb.get(i).pos.getPointer();
		}

		// convert the particle number to buffer size
		for (size_t i = 0; i < prc_sz_r.size(); i++)
		{
			prc_sz_r.get(i) = prc_sz_r.get(i) * (sizeof(prp_object) + sizeof(Point<dim, St> ));
		}

		// Send and receive the particles

		recv_mem_gm.clear();
		v_cl.sendrecvMultipleMessagesNBX(prc_sz_r.size(), (size_t *) prc_sz_r.getPointer(), (size_t *) prc_r.getPointer(), (void **) ptr.getPointer(), vector_dist_comm::message_alloc_map, this, NEED_ALL_SIZE);

		// Process the incoming particles

		process_received_map_list<prp_object, prp...>(v_pos,v_prp,out_part);

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
		// outgoing particles-id
		openfpm::vector<size_t> out_part;

		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// It contain the list of the processors this processor should to communicate with
		openfpm::vector<size_t> p_list;

		// map completely reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Contain the processor id of each particle (basically where they have to go)
		labelParticleProcessor<obp>(v_pos,opart, prc_sz, out_part);

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

		// Allocate the send buffers

		openfpm::vector<pos_prop> pb;

		// fill the send buffers
		fill_send_map_buf(v_pos,v_prp,prc_r, prc_sz_r, pb);

		// Create the set of pointers
		openfpm::vector<void *> ptr(prc_r.size());
		for (size_t i = 0; i < prc_r.size(); i++)
		{
			ptr.get(i) = pb.get(i).pos.getPointer();
		}

		// convert the particle number to buffer size
		for (size_t i = 0; i < prc_sz_r.size(); i++)
		{
			prc_sz_r.get(i) = prc_sz_r.get(i) * (sizeof(prop) + sizeof(Point<dim, St> ));
		}

		// Send and receive the particles

		recv_mem_gm.clear();
		v_cl.sendrecvMultipleMessagesNBX(prc_sz_r.size(), (size_t *) prc_sz_r.getPointer(), (size_t *) prc_r.getPointer(), (void **) ptr.getPointer(), vector_dist_comm::message_alloc_map, this, NEED_ALL_SIZE);

		// Process the incoming particles

		process_received_map(v_pos,v_prp,out_part);

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
		dec.swap(vc.dec);

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
		typedef openfpm::vector<prp_object, ExtPreAlloc<Memory>, typename memory_traits_lin<prp_object>::type, memory_traits_lin, openfpm::grow_policy_identity> send_vector;

		// Calculate memory and allocation for the send buffers

		// Total size
		size_t size_byte_prp = 0;
		size_t size_byte_pos = 0;

		calc_send_ghost_put_buf<prp_object>(v_pos,v_prp,size_byte_prp, size_byte_pos);

		// Create memory for the send buffer

		g_prp_mem.resize(size_byte_prp);

		// Create and fill send buffer for particle properties

		ExtPreAlloc<Memory> * prAlloc_prp = new ExtPreAlloc<Memory>(size_byte_prp, g_prp_mem);

		openfpm::vector<send_vector> g_send_prp;
		fill_send_ghost_put_prp_buf<send_vector, prp_object, prp...>(v_prp,g_send_prp, prAlloc_prp,g_m);

		// Send/receive the particle properties information
		v_cl.sendrecvMultipleMessagesNBX(prc_recv, g_send_prp, msg_alloc_ghost_put, this);
		process_received_put_ghost_prp<op,send_vector, prp_object, prp...>(v_prp,g_m);
	}
};


#endif /* SRC_VECTOR_VECTOR_DIST_COMM_HPP_ */
