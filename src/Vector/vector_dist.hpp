/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "HDF5_XdmfWriter/HDF5_XdmfWriter.hpp"
#include "VCluster.hpp"
#include "Space/Shape/Point.hpp"
#include "Vector/vector_dist_iterator.hpp"
#include "Space/Shape/Box.hpp"
#include "Vector/vector_dist_key.hpp"
#include "memory/PreAllocHeapMemory.hpp"
#include "memory/PtrMemory.hpp"
#include "NN/CellList/CellList.hpp"
#include "NN/CellList/CellListFast_hilb.hpp"
#include "util/common.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "CSVWriter/CSVWriter.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "Decomposition/common.hpp"
#include "Grid/grid_dist_id_iterator_dec.hpp"
#include "Grid/grid_key_dx_iterator_hilbert.hpp"
#include "Vector/vector_dist_ofb.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"

#define V_SUB_UNIT_FACTOR 64

#define NO_ID false
#define ID true

// Perform a ghost get or a ghost put
#define GET	1
#define PUT 2

#define NO_POSITION 1
#define WITH_POSITION 2

// Write the particles with ghost
#define NO_GHOST 0
#define WITH_GHOST 2

/*! \brief Distributed vector
 *
 * This class reppresent a distributed vector, the distribution of the structure
 * is based on the positional information of the elements the vector store
 *
 * ## Create a vector of random elements on each processor 2D
 * \snippet vector_dist_unit_test.hpp Create a vector of random elements on each processor 2D
 *
 * ## Create a vector of random elements on each processor 3D
 * \snippet vector_dist_unit_test.hpp Create a vector of random elements on each processor 3D
 *
 * ## Create a vector of elements distributed on a grid like way
 * \snippet vector_dist_unit_test.hpp Create a vector of elements distributed on a grid like way
 *
 * ## Redistribute the particles and sync the ghost properties
 * \snippet vector_dist_unit_test.hpp Redistribute the particles and sync the ghost properties
 *
 * \tparam dim Dimensionality of the space where the elements lives
 * \tparam St type of space float, double ...
 * \tparam prop properties the vector element store in OpenFPM data structure format
 * \tparam Decomposition Decomposition strategy to use CartDecomposition ...
 * \tparam Memory Memory pool where store the information HeapMemory ...
 *
 */

template<unsigned int dim, typename St, typename prop, typename Decomposition = CartDecomposition<dim,St>, typename Memory = HeapMemory>
class vector_dist
{
private:

	//! Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	size_t g_m = 0;

	//! Space Decomposition
	Decomposition dec;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	openfpm::vector<Point<dim, St>> v_pos;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	openfpm::vector<prop> v_prp;

	//! Virtual cluster
	Vcluster & v_cl;

	// definition of the send vector for position
	typedef openfpm::vector<Point<dim, St>, ExtPreAlloc<Memory>, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin , openfpm::grow_policy_identity> send_pos_vector;

	//////////////////////////////
	// COMMUNICATION variables
	//////////////////////////////

	//! It map the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! For each near processor, outgoing particle id and shift vector
	openfpm::vector<openfpm::vector<size_t>> opart;

	//! For each near processor, particle shift vector
	openfpm::vector<openfpm::vector<size_t>> oshift;

	//! For each adjacent processor the size of the ghost sending buffer
	openfpm::vector<size_t> ghost_prc_sz;

	//! Sending buffer for the ghost particles properties
	BHeapMemory g_prp_mem;

	//! Sending buffer for the ghost particles position
	BHeapMemory g_pos_mem;

	//! For each adjacent processor it store the size of the receiving message in byte
	openfpm::vector<size_t> recv_sz;

	//! For each adjacent processor it store the received message for ghost get
	openfpm::vector<BHeapMemory> recv_mem_gg;

	//! For each processor it store the received message for global map
	openfpm::vector<BHeapMemory> recv_mem_gm;

	/*! \brief It store for each processor the position and properties vector of the particles
	 *
	 *
	 */
	struct pos_prop
	{
		//! position vector
		openfpm::vector<Point<dim, St>, PreAllocHeapMemory<2>, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity> pos;
		//! properties vector
		openfpm::vector<prop, PreAllocHeapMemory<2>, typename memory_traits_lin<prop>::type, memory_traits_lin, openfpm::grow_policy_identity> prp;
	};

	template <typename sel_prop>
	struct pos_prop_sel
	{
		//! position vector
		openfpm::vector<Point<dim, St>, PreAllocHeapMemory<2>, typename memory_traits_lin<Point<dim, St>>::type, memory_traits_lin, openfpm::grow_policy_identity> pos;
		//! properties vector
		openfpm::vector<sel_prop, PreAllocHeapMemory<2>, typename memory_traits_lin<sel_prop>::type, memory_traits_lin, openfpm::grow_policy_identity> prp;
	};

	/*! \brief Label particles for mappings
	 *
	 * \param lbl_p Particle labeled
	 * \param prc_sz For each processor the number of particles to send
	 * \param opart id of the particles to send
	 *
	 */
	template<typename obp> void labelParticleProcessor(openfpm::vector<openfpm::vector<size_t>> & lbl_p, openfpm::vector<size_t> & prc_sz, openfpm::vector<size_t> & opart)
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
					prc_sz.get(p_id)++;lbl_p
					.get(p_id).add(key);
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

	/*! \brief Label the particles
	 *
	 * It count the number of particle to send to each processors and save its ids
	 *
	 * \param prc_sz For each processor the number of particles to send
	 * \param opart id if of the particles to send
	 * \param shift_id shift correction id
	 *
	 * \see nn_prcs::getShiftvectors()
	 *
	 */
	void labelParticlesGhost()
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
	 */
	void add_loc_particles_bc()
	{
		// get the shift vectors
		const openfpm::vector<Point<dim, St>> & shifts = dec.getShiftVectors();

		// this map is used to check if a combination is already present
		std::unordered_map<size_t, size_t> map_cmb;

		// Add local particles coming from periodic boundary, the only boxes that count are the one
		// touching the border, filter them

		// The boxes touching the border of the domain are divided in groups (first vector)
		// each group contain internal ghost coming from sub-domain of the same section
		openfpm::vector_std<openfpm::vector_std<Box<dim, St>>>box_f;
		openfpm::vector_std<comb<dim>> box_cmb;

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

		if (box_f.size() == 0)
			return;
		else
		{
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
							Point<dim, St> p = v_pos.get(key);
							// shift
							p -= shifts.get(box_cmb.get(i).lin());

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
	}

	/*! \brief This function fill the send buffer for the particle position after the particles has been label with labelParticles
	 *
	 * \param g_pos_send Send buffer to fill
	 * \param prAlloc_pos Memory object for the send buffer
	 *
	 */
	void fill_send_ghost_pos_buf(openfpm::vector<send_pos_vector> & g_pos_send, ExtPreAlloc<Memory> * prAlloc_pos)
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

	/*! \brief This function fill the send buffer for properties after the particles has been label with labelParticles
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param g_send_prp Send buffer to fill
	 * \param prAlloc_prp Memory object for the send buffer
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp> void fill_send_ghost_prp_buf(openfpm::vector<send_vector> & g_send_prp, ExtPreAlloc<Memory> * prAlloc_prp)
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
	 * \param prc_r List of processor rank involved in the send
	 * \param prc_r_sz For each processor in the list the size of the message to send
	 * \param pb send buffer
	 *
	 */
	void fill_send_map_buf(openfpm::vector<size_t> & prc_r, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<pos_prop> & pb)
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
	 * \param prc_r List of processor rank involved in the send
	 * \param prc_r_sz For each processor in the list the size of the message to send
	 * \param pb send buffer
	 *
	 */
	template<typename prp_object,int ... prp> void fill_send_map_buf_list(openfpm::vector<size_t> & prc_r, openfpm::vector<size_t> & prc_sz_r, openfpm::vector<pos_prop_sel<prp_object>> & pb)
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

	/*! \brief This function process the receiced data for the properties and populate the ghost
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp> void process_received_ghost_prp()
	{
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

			// resize with the number of elements
			v2.resize(n_ele);

			// Add the ghost particle
			v_prp.template add_prp<prp_object, PtrMemory, openfpm::grow_policy_identity, prp...>(v2);
		}
	}

	/*! \brief This function process the received data for the properties and populate the ghost
	 *
	 */
	void process_received_ghost_pos()
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
	 * \param list of the out-going particles
	 *
	 */
	void process_received_map(openfpm::vector<size_t> & out_part)
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
	 * \param list of the out-going particles
	 *
	 */
	template<typename prp_object , int ... prp> void process_received_map_list(openfpm::vector<size_t> & out_part)
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
	 * \param size_byte_prp total size for the property buffer
	 * \param size_byte_pos total size for the position buffer
	 *
	 */
	template<typename prp_object> void calc_send_ghost_buf(size_t & size_byte_prp, size_t & size_byte_pos)
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

	/*! \brief Calculate parameters for the cell list
	 *
	 * \param div Division array
	 * \param r_cut interation radius or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost. This parameter says how much must be enlarged
	 *
	 * \return the processor bounding box
	 */
	inline Box<dim, St> cl_param_calculate(size_t (&div)[dim], St r_cut, const Ghost<dim, St> & enlarge)

	{
		// calculate the parameters of the cell list

		// get the processor bounding box
		Box<dim, St> pbox = dec.getProcessorBounds();

		// extend by the ghost
		pbox.enlarge(enlarge);

		// Calculate the division array and the cell box
		for (size_t i = 0; i < dim; i++)
		{
			div[i] = static_cast<size_t>((pbox.getP2().get(i) - pbox.getP1().get(i)) / r_cut);
			div[i]++;
			pbox.setHigh(i,pbox.getLow(i) + div[i]*r_cut);
		}
		return pbox;
	}

	/*! \brief Initialize the structures
	 *
	 * \param np number of particles
	 *
	 */
	void init_structures(size_t np)
	{
		// convert to a local number of elements
		size_t p_np = np / v_cl.getProcessingUnits();

		// Get non divisible part
		size_t r = np % v_cl.getProcessingUnits();

		// Distribute the remain particles
		if (v_cl.getProcessUnitID() < r)
			p_np++;

		// resize the position vector
		v_pos.resize(p_np);

		// resize the properties vector
		v_prp.resize(p_np);

		g_m = p_np;
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

public:

	//! space type
	typedef St stype;

	//! dimensions of space
	static const unsigned int dims = dim;


	/*! \brief Constructor
	 *
	 * \param np number of elements
	 * \param box domain where the vector of elements live
	 * \param boundary conditions
	 * \param g Ghost margins
	 *
	 */
	vector_dist(const Decomposition & dec, size_t np) :
			dec(dec), v_cl(create_vcluster())
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		init_structures(np);
	}


	/*! \brief Constructor
	 *
	 * \param np number of elements
	 * \param box domain where the vector of elements live
	 * \param boundary conditions
	 * \param g Ghost margins
	 *
	 */
	vector_dist(size_t np, Box<dim, St> box, const size_t (&bc)[dim], const Ghost<dim, St> & g) :
			dec(create_vcluster()), v_cl(create_vcluster())
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		init_structures(np);
		init_decomposition(box,bc,g);
	}

	~vector_dist()
	{
#ifdef SE_CLASS2
		check_delete(this);
#endif
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

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local()
	{
		return g_m;
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPos(vect_dist_key_dx vec_key) -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getPos(vect_dist_key_dx vec_key) const -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
		return v_pos.template get<0>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) const -> const decltype(v_prp.template get<id>(vec_key.getKey()))
	{
		return v_prp.template get<id>(vec_key.getKey());
	}

	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * \tparam out of bound policy it specify what to do when the particles are detected out of bound
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particles
	 *
	 */
	template<typename obp = KillParticle> void map()
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
		labelParticleProcessor<obp>(opart, prc_sz, out_part);

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
		fill_send_map_buf(prc_r, prc_sz_r, pb);

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
		v_cl.sendrecvMultipleMessagesNBX(prc_sz_r.size(), (size_t *) prc_sz_r.getPointer(), (size_t *) prc_r.getPointer(), (void **) ptr.getPointer(), vector_dist::message_alloc_map, this, NEED_ALL_SIZE);

		// Process the incoming particles

		process_received_map(out_part);

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
	 */
	template<typename obp = KillParticle,unsigned int ... prp> void map_list()
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
		labelParticleProcessor<obp>(opart, prc_sz, out_part);

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
		fill_send_map_buf_list<prp_object,prp...>(prc_r, prc_sz_r, pb);

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
		v_cl.sendrecvMultipleMessagesNBX(prc_sz_r.size(), (size_t *) prc_sz_r.getPointer(), (size_t *) prc_r.getPointer(), (void **) ptr.getPointer(), vector_dist::message_alloc_map, this, NEED_ALL_SIZE);

		// Process the incoming particles

		process_received_map_list<prp_object, prp...>(out_part);

		// mark the ghost part

		g_m = v_pos.size();
	}

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam prp list of properties to get synchronize
	 * \param opt options WITH_POSITION, it send also the positional information of the particles
	 *
	 */
	template<int ... prp> void ghost_get(size_t opt = WITH_POSITION)
	{
		// Unload receive buffer
		for (size_t i = 0 ; i < recv_mem_gg.size() ; i++)
			recv_sz.get(i) = 0;

		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object, ExtPreAlloc<Memory>, typename memory_traits_lin<prp_object>::type, memory_traits_lin, openfpm::grow_policy_identity> send_vector;

		// reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Label all the particles
		labelParticlesGhost();

		// Calculate memory and allocation for the send buffers

		// Total size
		size_t size_byte_prp = 0;
		size_t size_byte_pos = 0;

		calc_send_ghost_buf<prp_object>(size_byte_prp, size_byte_pos);

		// Create memory for the send buffer

		g_prp_mem.resize(size_byte_prp);
		if (opt != NO_POSITION)
			g_pos_mem.resize(size_byte_pos);

		// Create and fill send buffer for particle properties

		ExtPreAlloc<Memory> * prAlloc_prp = new ExtPreAlloc<Memory>(size_byte_prp, g_prp_mem);

		openfpm::vector<send_vector> g_send_prp;
		fill_send_ghost_prp_buf<send_vector, prp_object, prp...>(g_send_prp, prAlloc_prp);

		// Create and fill the send buffer for the particle position

		ExtPreAlloc<Memory> * prAlloc_pos;
		openfpm::vector<send_pos_vector> g_pos_send;
		if (opt != NO_POSITION)
		{
			prAlloc_pos = new ExtPreAlloc<Memory>(size_byte_pos, g_pos_mem);
			fill_send_ghost_pos_buf(g_pos_send, prAlloc_pos);
		}

		// Create processor list
		openfpm::vector<size_t> prc;
		for (size_t i = 0; i < opart.size(); i++)
			prc.add(dec.IDtoProc(i));

		// Send/receive the particle properties information
		v_cl.sendrecvMultipleMessagesNBX(prc, g_send_prp, msg_alloc_ghost_get, this);
		process_received_ghost_prp<send_vector, prp_object, prp...>();

		if (opt != NO_POSITION)
		{
			// Send/receive the particle properties information
			v_cl.sendrecvMultipleMessagesNBX(prc, g_pos_send, msg_alloc_ghost_get, this);
			process_received_ghost_pos();
		}

		add_loc_particles_bc();
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
		vector_dist<dim, St, prop, Decomposition, Memory> * v = static_cast<vector_dist<dim, St, prop, Decomposition, Memory> *>(ptr);

		v->recv_sz.resize(v->dec.getNNProcessors());
		v->recv_mem_gg.resize(v->dec.getNNProcessors());

		// Get the local processor id
		size_t lc_id = v->dec.ProctoID(i);

		// resize the receive buffer
		v->recv_mem_gg.get(lc_id).resize(msg_i);
		v->recv_sz.get(lc_id) = msg_i;

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
		vector_dist<dim, St, prop, Decomposition, Memory> * vd = static_cast<vector_dist<dim, St, prop, Decomposition, Memory> *>(ptr);

		vd->recv_mem_gm.resize(vd->v_cl.getProcessingUnits());
		vd->recv_mem_gm.get(i).resize(msg_i);

		return vd->recv_mem_gm.get(i).getPointer();
	}

	/*! \brief Add local particle
	 *
	 * It add a local particle, with "local" we mean in this processor
	 * the particle can be also created out of the processor domain, in this
	 * case a call to map is required. Added particles are always created at the
	 * end and can be accessed with getLastPos and getLastProp
	 *
	 */
	void add()
	{
		v_prp.insert(g_m);
		v_pos.insert(g_m);

		g_m++;
	}

	/*! \brief Get the position of the last element
	 *
	 * \return the position of the element in space
	 *
	 */
	inline auto getLastPos() -> decltype(v_pos.template get<0>(0))
	{
		return v_pos.template get<0>(g_m - 1);
	}

	/*! \brief Get the property of the last element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> inline auto getLastProp() -> decltype(v_prp.template get<id>(0))
	{
		return v_prp.template get<id>(g_m - 1);
	}

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > CellL getCellList(St r_cut)
	{
		// Get ghost and anlarge by 1%
		Ghost<dim,St> g = dec.getGhost();
		g.magnify(1.013);

		return getCellList(r_cut, g);
	}

	/*! \brief Construct an hilbert cell list starting from the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 *
	 * \return the Cell list
	 *
	 */
	template<typename CellL = CellList_hilb<dim, St, FAST, shift<dim, St> > > CellL getCellList_hilb(St r_cut)
	{
		// Get ghost and anlarge by 1%
		Ghost<dim,St> g = dec.getGhost();
		g.magnify(1.013);

		return getCellList_hilb(r_cut, g);
	}

	/*! \brief Update a cell list using the stored particles
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param cell_list Cell list to update
	 *
	 */
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > void updateCellList(CellL & cell_list)
	{
		// Clear the cell list from the previous particles
		cell_list.clear();

		// for each particle add the particle to the cell list

		auto it = getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			cell_list.add(this->getPos(key), key.getKey());

			++it;
		}

		cell_list.set_gm(g_m);
	}

	/*! \brief Construct a cell list starting from the stored particles
	 *
	 * It differ from the get getCellList for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 *
	 */
	template<typename CellL = CellList<dim, St, FAST, shift<dim, St> > > CellL getCellList(St r_cut, const Ghost<dim, St> & enlarge)
	{
		CellL cell_list;

		// Division array
		size_t div[dim];

		// Processor bounding box
		auto pbox = cl_param_calculate(div, r_cut, enlarge);

		cell_list.Initialize(pbox, div);

		updateCellList(cell_list);

		return cell_list;
	}

	/*! \brief Construct an hilbert cell list starting from the stored particles
	 *
	 * It differ from the get getCellList for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param r_cut interation radius, or size of each cell
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 *
	 */
	template<typename CellL = CellList_hilb<dim, St, FAST, shift<dim, St> > > CellL getCellList_hilb(St r_cut, const Ghost<dim, St> & enlarge)
	{
		CellL cell_list;

		// Division array
		size_t div[dim];

		// Processor bounding box
		auto pbox = cl_param_calculate(div, r_cut, enlarge);

		cell_list.Initialize(pbox, div, g_m);

		updateCellList(cell_list);

		return cell_list;
	}

	/*! \brief for each particle get the verlet list
	 *
	 * \param verlet output verlet list for each particle
	 * \param r_cut cut-off radius
	 *
	 */
	void getVerlet(openfpm::vector<openfpm::vector<size_t>> & verlet, St r_cut)
	{
		// resize verlet to store the number of particles
		verlet.resize(size_local());

		// get the cell-list
		auto cl = getCellList(r_cut);

		// square of the cutting radius
		St r_cut2 = r_cut * r_cut;

		// iterate the particles
		auto it_p = this->getDomainIterator();
		while (it_p.isNext())
		{
			// key
			vect_dist_key_dx key = it_p.get();

			// Get the position of the particles
			Point<dim, St> p = this->getPos(key);

			// Clear the neighborhood of the particle
			verlet.get(key.getKey()).clear();

			// Get the neighborhood of the particle
			auto NN = cl.template getNNIterator<NO_CHECK>(cl.getCell(p));
			while (NN.isNext())
			{
				auto nnp = NN.get();

				// p != q
				if (nnp == key.getKey())
				{
					++NN;
					continue;
				}

				Point<dim, St> q = this->getPos(nnp);

				if (p.distance2(q) < r_cut2)
					verlet.get(key.getKey()).add(nnp);

				// Next particle
				++NN;
			}

			// next particle
			++it_p;
		}
	}

	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 * \tparam CellL CellList type to construct
	 *
	 * \param m an order of a hilbert curve
	 *
	 *
	 *
	 */
	template<typename CellL=CellList<dim,St,FAST,shift<dim,St> > > void reorder (int32_t m)
	{
		reorder(m,dec.getGhost());
	}


	/*! \brief Construct a cell list starting from the stored particles and reorder a vector according to the Hilberts curve
	 *
	 *
	 *It differs from the reorder(m) for an additional parameter, in case the
	 * domain + ghost is not big enough to contain additional padding particles, a Cell list
	 * with bigger space can be created
	 * (padding particles in general are particles added by the user out of the domains)
	 *
	 * \param m order of a curve
	 * \param enlarge In case of padding particles the cell list must be enlarged, like a ghost this parameter say how much must be enlarged
	 *
	 */
	template<typename CellL=CellList<dim,St,FAST,shift<dim,St> > > void reorder(int32_t m, const Ghost<dim,St> & enlarge)
	{
		// reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);


		CellL cell_list;

		// calculate the parameters of the cell list

		// get the processor bounding box
		Box<dim,St> pbox = dec.getProcessorBounds();
		// extend by the ghost
		pbox.enlarge(enlarge);

		size_t div[dim];

		// Calculate the division array and the cell box
		for (size_t i = 0 ; i < dim ; i++)
		{
			div[i] = 1 << m;
		}

		cell_list.Initialize(pbox,div);

		// for each particle add the particle to the cell list

		auto it = getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			cell_list.add(this->getPos(key),key.getKey());

			++it;
		}

		// Use cell_list to reorder v_pos

		//destination vector
		openfpm::vector<Point<dim,St>> v_pos_dest;
		openfpm::vector<prop> v_prp_dest;

		v_pos_dest.resize(v_pos.size());
		v_prp_dest.resize(v_prp.size());

		//hilberts curve iterator
		grid_key_dx_iterator_hilbert<dim> h_it(m);

		//Index for v_pos_dest
		size_t count = 0;

		grid_key_dx<dim> ksum;

		for (size_t i = 0; i < dim ; i++)
			ksum.set_d(i,cell_list.getPadding(i));

		while (h_it.isNext())
		{
		  auto key = h_it.get();
		  key += ksum;

		  size_t lin = cell_list.getGrid().LinId(key);

		  // for each particle in the Cell "lin"
		  for (size_t i = 0; i < cell_list.getNelements(lin); i++)
		  {
			  //reorder
			  auto v = cell_list.get(lin,i);
			  v_pos_dest.get(count) = v_pos.get(v);
			  v_prp_dest.get(count) = v_prp.get(v);

			  count++;
		  }
		  ++h_it;
		}

		v_pos.swap(v_pos_dest);
		v_prp.swap(v_prp_dest);
	}

	/*! \brief It return the number of particles contained by the previous processors
	 *
	 * \Warning It only work with the initial decomposition
	 *
	 * Given 1000 particles and 3 processors, you will get
	 *
	 * * Processor 0: 0
	 * * Processor 1: 334
	 * * Processor 2: 667
	 *
	 * \param np initial number of particles
	 *
	 */
	size_t init_size_accum(size_t np)
	{
		size_t accum = 0;

		// convert to a local number of elements
		size_t p_np = np / v_cl.getProcessingUnits();

		// Get non divisible part
		size_t r = np % v_cl.getProcessingUnits();

		accum = p_np * v_cl.getProcessUnitID();

		// Distribute the remain particles
		if (v_cl.getProcessUnitID() <= r)
			accum += v_cl.getProcessUnitID();
		else
			accum += r;

		return accum;
	}

	/*! \brief Get an iterator that traverse domain and ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getIterator()
	{
		return vector_dist_iterator(0, v_pos.size());
	}

	/*! /brief Get a grid Iterator
	 *
	 * Usefull function to place particles on a grid or grid-like (grid + noise)
	 *
	 * \return a Grid iterator
	 *
	 */
	inline grid_dist_id_iterator_dec<Decomposition> getGridIterator(const size_t (&sz)[dim])
	{
		size_t sz_g[dim];
		grid_key_dx<dim> start;
		grid_key_dx<dim> stop;
		for (size_t i = 0; i < dim; i++)
		{
			start.set_d(i, 0);
			if (dec.periodicity(i) == PERIODIC)
			{
				sz_g[i] = sz[i];
				stop.set_d(i, sz_g[i] - 2);
			}
			else
			{
				sz_g[i] = sz[i];
				stop.set_d(i, sz_g[i] - 1);
			}
		}

		grid_dist_id_iterator_dec<Decomposition> it_dec(dec, sz_g, start, stop);
		return it_dec;
	}

	void calculateComputationCosts()
	{

	}

	/*! \brief Get the iterator across the position of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getGhostIterator() const
	{
		return vector_dist_iterator(g_m, v_pos.size());
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator getDomainIterator() const
	{
		return vector_dist_iterator(0, g_m);
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	Decomposition & getDecomposition()
	{
		return dec;
	}

	/*! \brief Remove a set of elements from the distributed vector
	 *
	 * \warning keys must be sorted
	 *
	 * \param keys vector of elements to eliminate
	 * \param start from where to eliminate
	 *
	 */
	void remove(openfpm::vector<size_t> & keys, size_t start = 0)
	{
		v_pos.remove(keys, start);
		v_prp.remove(keys, start);

		g_m -= keys.size();
	}

	/*! \brief Remove one element from the distributed vector
	 *
	 * \warning keys must be sorted
	 *
	 * \param keys vector of elements to eliminate
	 * \param start from where to eliminate
	 *
	 */
	void remove(size_t key)
	{
		v_pos.remove(key);
		v_prp.remove(key);

		g_m--;
	}

	inline void addComputationCosts()
	{
		CellDecomposer_sm<dim, St> cdsm;

		cdsm.setDimensions(dec.getDomain(), dec.getGrid().getSize(), 0);

		for (size_t i = 0; i < dec.getNSubSubDomains(); i++)
		{
			dec.setSubSubDomainComputationCost(i, 1);
		}

		auto it = getDomainIterator();

		while (it.isNext())
		{
			size_t v = cdsm.getCell(this->getPos(it.get()));

			dec.addComputationCost(v, 1);

			++it;
		}

	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param opt VTK_WRITER or CSV_WRITER
	 *
	 * \return true if the file has been written without error
	 *
	 */
	inline bool write(std::string out, int opt = NO_GHOST | VTK_WRITER )
	{

		if ((opt & 0xFFFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim,St>>, openfpm::vector<prop> > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output,v_pos,v_prp);
		}
		else if ((opt & 0xFFFF0000) == VTK_WRITER)
		{
			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,St>>, openfpm::vector<prop>>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output);
		}

		return false;
	}

	/*! \brief Delete the particles on the ghost
	 *
	 *
	 */
	void deleteGhost()
	{
		v_pos.resize(g_m);
		v_prp.resize(g_m);
	}

	/*! \brief Resize the vector (locally)
	 *
	 * \warning It automaticallt delete the ghosts
	 *
	 * \param rs
	 *
	 */
	void resize(size_t rs)
	{
		deleteGhost();

		v_pos.resize(rs);
		v_prp.resize(rs);

		g_m = rs;
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param opt NO_GHOST or WITH_GHOST
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write(std::string out, size_t iteration, int opt = NO_GHOST)
	{
		if ((opt & 0xFFFF0000) == CSV_WRITER)
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<Point<dim, St>>, openfpm::vector<prop> > csv_writer;

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".csv"));

			// Write the CSV
			return csv_writer.write(output, v_pos, v_prp);
		}
		else
		{
			// VTKWriter for a set of points
			VTKWriter<boost::mpl::pair<openfpm::vector<Point<dim,St>>, openfpm::vector<prop>>, VECTOR_POINTS> vtk_writer;
			vtk_writer.add(v_pos,v_prp,g_m);

			std::string output = std::to_string(out + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(iteration) + std::to_string(".vtk"));

			// Write the VTK file
			return vtk_writer.write(output);
		}
	}

	/* \brief It return the id of structure in the allocation list
	 *
	 * \see print_alloc and SE_CLASS2
	 *
	 */
	long int who()
	{
#ifdef SE_CLASS2
		return check_whoami(this,8);
#else
		return -1;
#endif
	}

	/*! \brief Get the Virtual Cluster machine
	 *
	 * \return the Virtual cluster machine
	 *
	 */

	Vcluster & getVC()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return v_cl;
	}
};


#endif /* VECTOR_HPP_ */
