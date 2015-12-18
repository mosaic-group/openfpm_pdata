/*
 * Vector.hpp
 *
 *  Created on: Mar 5, 2015
 *      Author: Pietro Incardona
 */

#ifndef VECTOR_HPP_
#define VECTOR_HPP_

#include "VCluster.hpp"
#include "Space/Shape/Point.hpp"
#include "Vector/vector_dist_iterator.hpp"
#include "Space/Shape/Box.hpp"
#include "Vector/vector_dist_key.hpp"
#include "memory/PreAllocHeapMemory.hpp"
#include "memory/PtrMemory.hpp"
#include "NN/CellList/CellList.hpp"
#include "util/common.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "CSVWriter.hpp"
#include "Decomposition/common.hpp"

#define V_SUB_UNIT_FACTOR 64

#define NO_ID false
#define ID true

// Perform a ghost get or a ghost put
#define GET	1
#define PUT 2

#define INTERNAL 0

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

template<unsigned int dim, typename St, typename prop, typename Decomposition , typename Memory=HeapMemory>
class vector_dist
{
private:

	//! Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	size_t g_m = 0;

	//! Space Decomposition
	Decomposition dec;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	Vcluster_object_array<openfpm::vector<Point<dim,St>>> v_pos;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	Vcluster_object_array<openfpm::vector<prop>> v_prp;

	//! Virtual cluster
	Vcluster & v_cl;

public:

	/*! \brief Constructor
	 *
	 * \param np number of elements
	 * \param box domain where the vector of elements live
	 * \param g Ghost margins
	 *
	 */
	vector_dist(size_t np, Box<dim,St> box, Ghost<dim,St> g = Ghost<dim,St>())
	:dec(*global_v_cluster),v_cl(*global_v_cluster)
	{
#ifdef SE_CLASS2
		check_new(this,8,VECTOR_DIST_EVENT,4);
#endif

		// Allocate unassigned particles vectors
		v_pos = v_cl.template allocate<openfpm::vector<Point<dim,St>>>(1);
		v_prp = v_cl.template allocate<openfpm::vector<prop>>(1);

		// convert to a local number of elements
		size_t p_np = np / v_cl.getProcessingUnits();

		// Get non divisible part
		size_t r = np % v_cl.getProcessingUnits();

		// Distribute the remain particles
		if (v_cl.getProcessUnitID() < r)
			p_np++;

		// resize the position vector
		v_pos.get(0).resize(p_np);

		// resize the properties vector
		v_prp.get(0).resize(p_np);

		// Create a valid decomposition of the space
		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * getDefaultNsubsub();

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/dim));}

		// boundary conditions
		size_t bc[dim];
		for (size_t i = 0 ; i < dim ; i++)
			bc[i] = NON_PERIODIC;

		// Create the sub-domains
		dec.setParameters(div,box,bc,g);

		Point<dim,St> p;
		p.zero();
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
		return  V_SUB_UNIT_FACTOR;
	}

	/*! \brief return the local size of the vector
	 *
	 * \return local size
	 *
	 */
	size_t size_local()
	{
		return v_pos.get(0).size();
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
	template<unsigned int id> inline auto getPos(vect_dist_key_dx vec_key) -> decltype(v_pos.get(vec_key.getSub()).template get<id>(vec_key.getKey()))
	{
		return v_pos.get(vec_key.getSub()).template get<id>(vec_key.getKey());
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
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(v_prp.get(vec_key.getSub()).template get<id>(vec_key.getKey()))
	{
		return v_prp.get(vec_key.getSub()).template get<id>(vec_key.getKey());
	}

	/*! \brief It store for each processor the position and properties vector of the particles
	 *
	 *
	 */
	struct pos_prop
	{
		//! position vector
		openfpm::vector<Point<dim,St>,PreAllocHeapMemory<2>,openfpm::grow_policy_identity> pos;
		//! properties vector
		openfpm::vector<prop,PreAllocHeapMemory<2>,openfpm::grow_policy_identity> prp;
	};

	//! It map the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	/*! \brief It move all the particles that does not belong to the local processor to the respective processor
	 *
	 * In general this function is called after moving the particles to move the
	 * elements out the local processor. Or just after initialization if each processor
	 * contain non local particle
	 *
	 */
	void map()
	{
		// outgoing particles-id
		openfpm::vector<size_t> opart;

		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// Unassigned particle vector, is always the last vector
		size_t up_v = v_pos.size()-1;

		// Contain the map of the processors, this processors should communicate with
		openfpm::vector<size_t> p_map(v_cl.getProcessingUnits());

		// Contain the processor id of each particle (basically where they have to go)
		openfpm::vector<size_t> lbl_p(v_pos.get(up_v).size());

		// It contain the list of the processors this processor should to communicate with
		openfpm::vector<size_t> p_list;

		auto it = v_pos.get(up_v).getIterator();

		// Label all the particles with the processor id where they should go
		while (it.isNext())
		{
			auto key = it.get();

			size_t p_id = dec.processorID(v_pos.get(up_v).get(key));

			lbl_p.get(key) = p_id;

			// It has to communicate
			if (p_id != v_cl.getProcessUnitID())
			{
				p_map.get(p_id) = 1;
				prc_sz.get(p_id)++;

				opart.add(key);
			}

			// Add processors and add size

			++it;
		}

		// resize the map
		p_map_req.resize(v_cl.getProcessingUnits());

		// Create the sz and prc buffer

		openfpm::vector<size_t> prc_sz_r;
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0 ; i < v_cl.getProcessingUnits() ; i++)
		{
			if (p_map.get(i) == 1)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(prc_sz.get(i));
			}
		}

		// Allocate all the buffers

		openfpm::vector<pos_prop> pb(prc_r.size());

		for (size_t i = 0 ;  i < prc_r.size() ; i++)
		{
			// Create the size required to store the particles position and properties to communicate
			size_t s1 = openfpm::vector<Point<dim,St>,HeapMemory,openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i),0);
			size_t s2 = openfpm::vector<prop,HeapMemory,openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i),0);

			// Preallocate the memory
			size_t sz[2] = {s1,s2};
			PreAllocHeapMemory<2> * mem = new PreAllocHeapMemory<2>(sz);

			// Set the memory allocator
			pb.get(i).pos.setMemory(*mem);
			pb.get(i).prp.setMemory(*mem);

			// set the size and allocate, using mem warant that pos and prp is contiguous
			pb.get(i).pos.resize(prc_sz_r.get(i));
			pb.get(i).prp.resize(prc_sz_r.get(i));
		}


		// Run through all the particles and fill pb, the sending buffer

		openfpm::vector<size_t> prc_cnt(prc_r.size());
		prc_cnt.fill(0);

		it = lbl_p.getIterator();

		while (it.isNext())
		{
			auto key = it.get();

			size_t lbl = lbl_p.get(key);
			if (lbl == v_cl.getProcessUnitID())
			{
				++it;
				continue;
			}

			lbl = p_map_req.get(lbl);

			pb.get(lbl).pos.set(prc_cnt.get(lbl),v_pos.get(up_v).get(key));
			pb.get(lbl).prp.set(prc_cnt.get(lbl),v_prp.get(up_v).get(key));

			prc_cnt.get(lbl)++;

			// Add processors and add size
			++it;
		}

		// Create the set of pointers
		openfpm::vector<void *> ptr(prc_r.size());
		for (size_t i = 0 ; i < prc_r.size() ; i++)
		{
			ptr.get(i) = pb.get(i).pos.getPointer();
		}

		// convert the particle number to buffer size
		for (size_t i = 0 ; i < prc_sz_r.size() ; i++)
		{
			prc_sz_r.get(i) = prc_sz_r.get(i)*(sizeof(prop) + sizeof(Point<dim,St>));
		}

		// Send and receive the particles

		recv_cnt = 0;
		v_cl.sendrecvMultipleMessagesPCX(prc_sz_r.size(),&p_map.get(0), (size_t *)prc_sz_r.getPointer(), (size_t *)prc_r.getPointer() , (void **)ptr.getPointer() , vector_dist::message_alloc_map, this ,NEED_ALL_SIZE);

		// overwrite the outcoming particle with the incoming particle and resize the vectors

		size_t total_element = 0;
		size_t o_p_id = 0;

		for (size_t i = 0 ; i < v_proc.size() ; i++)
		{
			// Get the number of elements

			size_t n_ele = v_proc.get(i) / (sizeof(Point<dim,St>) + sizeof(prop));

			// Pointer of the received positions for each near processor
			void * ptr_pos = ((unsigned char *)hp_recv.getPointer()) + (total_element * (sizeof(Point<dim,St>) + sizeof(prop)));
			// Pointer of the received properties for each near processor
			void * ptr_prp = ((unsigned char *)hp_recv.getPointer()) + (total_element * (sizeof(Point<dim,St>) + sizeof(prop))) + n_ele * sizeof(Point<dim,St>);

			PtrMemory * ptr1 = new PtrMemory(ptr_pos,n_ele * sizeof(Point<dim,St>));
			PtrMemory * ptr2 = new PtrMemory(ptr_prp,n_ele * sizeof(prop));

			// create vector representation to a piece of memory already allocated

			openfpm::vector<Point<dim,St>,PtrMemory,openfpm::grow_policy_identity> vpos;
			openfpm::vector<prop,PtrMemory,openfpm::grow_policy_identity> vprp;

			vpos.setMemory(*ptr1);
			vprp.setMemory(*ptr2);

			vpos.resize(n_ele);
			vprp.resize(n_ele);

			// Add the received particles to v_pos and v_prp

			size_t j = 0;
			for ( ; j < vpos.size() && o_p_id < opart.size() ; j++, o_p_id++)
			{
				v_pos.get(0).set(opart.get(o_p_id),vpos.get(j));
				v_prp.get(0).set(opart.get(o_p_id),vprp.get(j));
			}

			for ( ; j < vpos.size(); j++)
			{
				v_pos.get(0).add();
				v_pos.get(0).set(v_pos.get(0).size()-1,vpos.get(j));
				v_prp.get(0).add();
				v_prp.get(0).set(v_prp.get(0).size()-1,vprp.get(j));
			}

			// increment the total number of element counter
			total_element += n_ele;
		}

		// remove the hole (out-going particles) in the vector

		v_pos.get(0).remove(opart,o_p_id);
		v_prp.get(0).remove(opart,o_p_id);
	}

	//! For each adjacent processor outgoing particles-ids
	openfpm::vector<openfpm::vector<size_t>> opart;

	//! For each adjacent processor the size of the ghost sending buffer
	openfpm::vector<size_t> ghost_prc_sz;

	//! Sending buffer for the ghost particles properties
	Memory g_prp_mem;

	//! Sending buffer for the ghost particles position
	Memory g_pos_mem;

	/*! \brief It synchronize the properties and position of the ghost particles
	 *
	 * \tparam prp list of properties to get synchronize
	 * \param opt options WITH_POSITION, it send also the positional information of the particles
	 *
	 */
	template<int... prp> void ghost_get(size_t opt = WITH_POSITION)
	{
		// Sending property object
		typedef object<typename object_creator<typename prop::type,prp...>::type> prp_object;

		// send vector for each processor
		typedef  openfpm::vector<prp_object,ExtPreAlloc<Memory>> send_vector;

		// Buffer that contain the number of elements to send for each processor
		ghost_prc_sz.clear();
		ghost_prc_sz.resize(dec.getNNProcessors());
		// Buffer that contain for each processor the id of the particle to send
		opart.clear();
		opart.resize(dec.getNNProcessors());

		// Label the internal (assigned) particles
		auto it = v_pos.get(INTERNAL).getIterator();

		// Label all the particles with the processor id, where they should go
		while (it.isNext())
		{
			auto key = it.get();

			const openfpm::vector<size_t> & vp_id = dec.template ghost_processorID<typename Decomposition::lc_processor_id>(v_pos.get(INTERNAL).get(key),UNIQUE);

			for (size_t i = 0 ; i < vp_id.size() ; i++)
			{
				// processor id
				size_t p_id = vp_id.get(i);

				// add particle to communicate
				ghost_prc_sz.get(p_id)++;

				opart.get(p_id).add(key);
			}

			++it;
		}

		// Send buffer size in byte ( one buffer for all processors )
		size_t size_byte_prp = 0;
		size_t size_byte_pos = 0;

		// sequence of pre-allocation pattern for property and position send buffer
		std::vector<size_t> pap_prp;
		std::vector<size_t> pap_pos;

		// Calculate the total size required for the sending buffer
		for ( size_t i = 0 ; i < ghost_prc_sz.size() ; i++ )
		{
			size_t alloc_ele = openfpm::vector<prp_object>::calculateMem(ghost_prc_sz.get(i),0);
			pap_prp.push_back(alloc_ele);
			size_byte_prp += alloc_ele;

			alloc_ele = openfpm::vector<Point<dim,St>>::calculateMem(ghost_prc_sz.get(i),0);
			pap_pos.push_back(alloc_ele);
			size_byte_pos += alloc_ele;
		}

		// resize the property buffer memory
		g_prp_mem.resize(size_byte_prp);
		// resize the position buffer memory
		if (opt != NO_POSITION) g_pos_mem.resize(size_byte_pos);

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> * prAlloc_prp = new ExtPreAlloc<Memory>(pap_prp,g_prp_mem);

		ExtPreAlloc<Memory> * prAlloc_pos;
		// Create an object of preallocated memory for position
		if (opt != NO_POSITION) prAlloc_pos = new ExtPreAlloc<Memory>(pap_pos,g_pos_mem);

		// create a vector of send vector (ExtPreAlloc warrant that all the created vector are contiguous)
		openfpm::vector<send_vector> g_send_prp;

		// create a number of send buffers equal to the near processors
		g_send_prp.resize(ghost_prc_sz.size());
		for (size_t i = 0 ; i < g_send_prp.size() ; i++)
		{
			// set the preallocated memory to ensure contiguity
			g_send_prp.get(i).setMemory(*prAlloc_prp);

			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(ghost_prc_sz.get(i));
		}

		// Fill the send buffer
		for ( size_t i = 0 ; i < opart.size() ; i++ )
		{
			for (size_t j = 0 ; j < opart.get(i).size() ; j++)
			{
				// source object type
				typedef encapc<1,prop,typename openfpm::vector<prop>::memory_conf> encap_src;
				// destination object type
				typedef encapc<1,prp_object,typename openfpm::vector<prp_object>::memory_conf> encap_dst;

				// Copy only the selected properties
				object_si_d<encap_src,encap_dst,OBJ_ENCAP,prp...>(v_prp.get(INTERNAL).get(opart.get(i).get(j)),g_send_prp.get(i).get(j));
			}
		}

		// Create the buffer for particle position

		// definition of the send vector for position for each processor
		typedef  openfpm::vector<Point<dim,St>,ExtPreAlloc<Memory>> send_pos_vector;

		openfpm::vector<send_pos_vector> g_pos_send;
		if (opt != NO_POSITION)
		{
			// create a number of send buffers equal to the near processors
			g_pos_send.resize(ghost_prc_sz.size());
			for (size_t i = 0 ; i < g_pos_send.size() ; i++)
			{
				// set the preallocated memory to ensure contiguity
				g_pos_send.get(i).setMemory(*prAlloc_pos);

				// resize the sending vector (No allocation is produced)
				g_pos_send.get(i).resize(ghost_prc_sz.get(i));
			}

			// Fill the send buffer
			for ( size_t i = 0 ; i < opart.size() ; i++ )
			{
				for (size_t j = 0 ; j < opart.get(i).size() ; j++)
				{
					g_pos_send.get(i).set(j,v_pos.get(INTERNAL).get(opart.get(i).get(j)));
				}
			}
		}

		// Create processor buffer pattern

		openfpm::vector<size_t> prc;
		for (size_t i = 0 ; i < opart.size() ; i++)
		{
			prc.add(dec.IDtoProc(i));
		}

		// Send receive the particles properties information
		v_cl.sendrecvMultipleMessagesNBX(prc,g_send_prp,msg_alloc_ghost_get,this);

		// Mark the ghost part
		g_m = v_prp.get(INTERNAL).size();

		// Process the received data (recv_mem_gg != 0 if you have data)
		for (size_t i = 0 ; i < dec.getNNProcessors() && recv_mem_gg.size() != 0 ; i++)
		{
			// calculate the number of received elements
			size_t n_ele = recv_sz.get(i) / sizeof(prp_object);

			// add the received particles to the vector
			PtrMemory * ptr1 = new PtrMemory(recv_mem_gg.get(i).getPointer(),recv_sz.get(i));

			// create vector representation to a piece of memory already allocated
			openfpm::vector<prp_object,PtrMemory,openfpm::grow_policy_identity> v2;

			v2.setMemory(*ptr1);

			// resize with the number of elements
			v2.resize(n_ele);

			// Add the ghost particle
			v_prp.get(INTERNAL).template add_prp<prp_object,PtrMemory,openfpm::grow_policy_identity,prp...>(v2);
		}

		if (opt != NO_POSITION)
		{
			// Send receive the particles properties information
			v_cl.sendrecvMultipleMessagesNBX(prc,g_pos_send,msg_alloc_ghost_get,this);

			// Process the received data (recv_mem_gg != 0 if you have data)
			for (size_t i = 0 ; i < dec.getNNProcessors() && recv_mem_gg.size() != 0 ; i++)
			{
				// calculate the number of received elements
				size_t n_ele = recv_sz.get(i) / sizeof(Point<dim,St>);

				// add the received particles to the vector
				PtrMemory * ptr1 = new PtrMemory(recv_mem_gg.get(i).getPointer(),recv_sz.get(i));

				// create vector representation to a piece of memory already allocated

				openfpm::vector<Point<dim,St>,PtrMemory,openfpm::grow_policy_identity> v2;

				v2.setMemory(*ptr1);

				// resize with the number of elements
				v2.resize(n_ele);

				// Add the ghost particle
				v_pos.get(INTERNAL).template add<PtrMemory,openfpm::grow_policy_identity>(v2);
			}
		}
	}

	//! For each adjacent processor it store the size of the receiving message in byte
	openfpm::vector<size_t> recv_sz;

	//! For each adjacent processot it store the received message
	openfpm::vector<HeapMemory> recv_mem_gg;

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
	static void * msg_alloc_ghost_get(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		vector_dist<dim,St,prop,Decomposition,Memory> * v = static_cast<vector_dist<dim,St,prop,Decomposition,Memory> *>(ptr);

		v->recv_sz.resize(v->dec.getNNProcessors());
		v->recv_mem_gg.resize(v->dec.getNNProcessors());

		// Get the local processor id
		size_t lc_id = v->dec.ProctoID(i);

		// resize the receive buffer
		v->recv_mem_gg.get(lc_id).resize(msg_i);
		v->recv_sz.get(lc_id) = msg_i;

		return v->recv_mem_gg.get(lc_id).getPointer();
	}

	//! Receive buffer for global communication
	HeapMemory hp_recv;

	//! For each message contain the processor from which processor come from
	openfpm::vector<size_t> v_proc;

	//! Total size of the received buffer
	size_t recv_cnt;

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
	static void * message_alloc_map(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		// cast the pointer
		vector_dist<dim,St,prop,Decomposition,Memory> * vd = static_cast<vector_dist<dim,St,prop,Decomposition,Memory> *>(ptr);

		// Resize the receive buffer, and the size of each message buffer
		vd->hp_recv.resize(total_msg);
		vd->v_proc.resize(total_p);

		// Return the receive pointer
		void * recv_ptr = (unsigned char *)vd->hp_recv.getPointer() + vd->recv_cnt;

		// increment the receive pointer
		vd->recv_cnt += msg_i;

		// Save the processor message size
		vd->v_proc.get(ri) = msg_i;

		return recv_ptr;
	}


	/*! \brief Get the iterator across the position of the particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<Point<dim,St>>> getIterator()
	{
		return vector_dist_iterator<openfpm::vector<Point<dim,St>>>(v_pos);
	}

	/*! \brief Get the iterator across the position of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<Point<dim,St>>> getGhostIterator()
	{
		return vector_dist_iterator<openfpm::vector<Point<dim,St>>>(v_pos,g_m);
	}

	/*! \brief Get the iterator across the properties of the particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<prop>> getPropIterator()
	{
		return vector_dist_iterator<openfpm::vector<prop>>(v_prp);
	}

	/*! \brief Get the iterator across the properties of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<prop>> getGhostPropIterator()
	{
		return vector_dist_iterator<openfpm::vector<prop>>(v_prp,g_m);
	}

	/*! \brief Get the decomposition
	 *
	 * \return
	 *
	 */
	const Decomposition & getDecomposition()
	{
		return dec;
	}

	/*! \brief Output particle position and properties
	 *
	 * \param out output
	 * \param opt NO_GHOST or WITH_GHOST
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write(std::string out, int opt = NO_GHOST)
	{
		// CSVWriter test
		CSVWriter<openfpm::vector<Point<dim,St>>, openfpm::vector<prop> > csv_writer;

		std::string output = std::to_string(out + std::to_string(v_cl.getProcessUnitID()) + std::to_string(".csv"));

		// Write the CSV
		return csv_writer.write(output,v_pos.get(INTERNAL),v_prp.get(INTERNAL));
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
};


#endif /* VECTOR_HPP_ */
