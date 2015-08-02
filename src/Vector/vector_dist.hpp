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
 */

template<typename point, typename prop, typename Box, typename Decomposition , typename Memory=HeapMemory, bool with_id=false>
class vector_dist
{
private:

	// Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	size_t g_m = 0;

	// indicate from where the ghost particle start in the vector
	size_t ghost_pointer;

	//! Space Decomposition
	Decomposition dec;

	// Particle position vector for each sub-domain the last one is the unassigned particles vector
	Vcluster_object_array<openfpm::vector<point>> v_pos;

	// Particle properties vector for each sub-domain the last one is the unassigned particles vector
	Vcluster_object_array<openfpm::vector<prop>> v_prp;

	// Virtual cluster
	Vcluster & v_cl;

	// Geometrical cell list
	CellList<point::dims,typename point::coord_type,FAST> geo_cell;

	// Label particles


public:

	/*! \brief Constructor
	 *
	 * \param Global number of elements
	 *
	 */
	vector_dist(size_t np, Box box, Ghost<point::dims,typename point::coord_type> g = Ghost<point::dims,typename point::coord_type>())
	:dec(Decomposition(*global_v_cluster)),v_cl(*global_v_cluster)
	{
		// Allocate unassigned particles vectors
		v_pos = v_cl.template allocate<openfpm::vector<point>>(1);
		v_prp = v_cl.template allocate<openfpm::vector<prop>>(1);

		// convert to a local number of elements
		np /= v_cl.getProcessingUnits();

		// resize the position vector
		v_pos.get(0).resize(np);

		// resize the properties vector
		v_prp.get(0).resize(np);

		// Create a valid decomposition of the space
		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * SUB_UNIT_FACTOR;

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[point::dims];
		for (size_t i = 0 ; i < point::dims ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/point::dims));}

		// Create the sub-domains
		dec.setParameters(div,box,g);

		// Get the bounding box containing the processor domain
		const ::Box<point::dims,typename point::coord_type> & bbound = dec.getProcessorBounds();

		const ::Box<point::dims,typename point::coord_type> & smallest_unit = dec.getSmallestSubdivision();

		// convert spacing divisions
		size_t n_g[point::dims];

		for (size_t i = 0 ; i < point::dims ; i++)
			n_g[i] = (bbound.getHigh(i) - bbound.getLow(i)) / smallest_unit.getHigh(i);

		point p;
		p.zero();

		// Initialize the geo cell list
		geo_cell.Initialize(box,n_g,p,8);
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

	/*! \brief Get position of an object
	 *
	 * \param vec_key vector element
	 *
	 */
	template<unsigned int id> inline auto getPos(vect_dist_key_dx vec_key) -> decltype(v_pos.get(vec_key.getSub()).template get<id>(vec_key.getKey()))
	{
		return v_pos.get(vec_key.getSub()).template get<id>(vec_key.getKey());
	}

	/*! \brief Get the property of the object
	 *
	 * \param vec_key vector element
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
		openfpm::vector<point,openfpm::device_cpu<point>,PreAllocHeapMemory<2>,openfpm::grow_policy_identity> pos;
		//! properties vector
		openfpm::vector<prop,openfpm::device_cpu<prop>,PreAllocHeapMemory<2>,openfpm::grow_policy_identity> prp;
	};

	/*! \brief set the ghost
	 *
	 *  \param g ghost
	 *
	 */
	void setGhost()
	{
		dec.calculateGhostBoxes();
	}

	//! It map the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	/*! \brief It communicate the particle to the respective processor
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
			size_t s1 = openfpm::vector<point,openfpm::device_cpu<point>,HeapMemory,openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i),0);
			size_t s2 = openfpm::vector<prop,openfpm::device_cpu<prop>,HeapMemory,openfpm::grow_policy_identity>::calculateMem(prc_sz_r.get(i),0);

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
			prc_sz_r.get(i) = prc_sz_r.get(i)*(sizeof(prop) + sizeof(point));
		}

		// Send and receive the particles

		recv_cnt = 0;
		v_cl.sendrecvMultipleMessagesPCX(prc_sz_r.size(),&p_map.get(0), &prc_sz_r.get(0), &prc_r.get(0) , &ptr.get(0) , vector_dist::message_alloc_map, this ,NEED_ALL_SIZE);

		// overwrite the outcoming particle with the incoming particle and resize the vectors

		size_t total_element = 0;
		size_t o_p_id = 0;

		for (size_t i = 0 ; i < v_proc.size() ; i++)
		{
			// Get the number of elements

			size_t n_ele = v_proc.get(i) / (sizeof(point) + sizeof(prop));

			// Pointer of the received positions for each near processor
			void * ptr_pos = ((unsigned char *)hp_recv.getPointer()) + (total_element * (sizeof(point) + sizeof(prop)));
			// Pointer of the received properties for each near processor
			void * ptr_prp = ((unsigned char *)hp_recv.getPointer()) + (total_element * (sizeof(point) + sizeof(prop))) + n_ele * sizeof(point);

			PtrMemory * ptr1 = new PtrMemory(ptr_pos,n_ele * sizeof(point));
			PtrMemory * ptr2 = new PtrMemory(ptr_prp,n_ele * sizeof(prop));

			// create vector representation to a piece of memory already allocated

			openfpm::vector<point,openfpm::device_cpu<point>,PtrMemory,openfpm::grow_policy_identity> vpos;
			openfpm::vector<prop,openfpm::device_cpu<prop>,PtrMemory,openfpm::grow_policy_identity> vprp;

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

	// outgoing particles-id
	openfpm::vector<openfpm::vector<size_t>> opart;

	// Each entry contain the size of the ghost sending buffer
	openfpm::vector<size_t> ghost_prc_sz;

	// ghost particle labels
	openfpm::vector<size_t> ghost_lbl_p;

	// Memory for the ghost sending buffer
	Memory g_prp_mem;

	// Memory for the ghost position sending buffer
	Memory g_pos_mem;

	/*! \brief It synchronize getting the ghost particles
	 *
	 * \prp Properties to get
	 * \opt options WITH_POSITION, it send also the positional information of the particles
	 *
	 */
	template<int... prp> void ghost_get(size_t opt = WITH_POSITION)
	{
		// Sending property object
		typedef object<typename object_creator<typename prop::type,prp...>::type> prp_object;

		// send vector for each processor
		typedef  openfpm::vector<prp_object,openfpm::device_cpu<prp_object>,ExtPreAlloc<Memory>> send_vector;

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

			alloc_ele = openfpm::vector<point>::calculateMem(ghost_prc_sz.get(i),0);
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
				object_si_d<encap_src,encap_dst,ENCAP,prp...>(v_prp.get(INTERNAL).get(opart.get(i).get(j)),g_send_prp.get(i).get(j));
			}
		}

		// Create the buffer for particle position

		// definition of the send vector for position for each processor
		typedef  openfpm::vector<point,openfpm::device_cpu<point>,ExtPreAlloc<Memory>> send_pos_vector;

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
			openfpm::vector<prp_object,openfpm::device_cpu<prp_object>,PtrMemory,openfpm::grow_policy_identity> v2;

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
				size_t n_ele = recv_sz.get(i) / sizeof(point);

				// add the received particles to the vector
				PtrMemory * ptr1 = new PtrMemory(recv_mem_gg.get(i).getPointer(),recv_sz.get(i));

				// create vector representation to a piece of memory already allocated

				openfpm::vector<point,openfpm::device_cpu<point>,PtrMemory,openfpm::grow_policy_identity> v2;

				v2.setMemory(*ptr1);

				// resize with the number of elements
				v2.resize(n_ele);

				// Add the ghost particle
				v_pos.get(INTERNAL).template add<PtrMemory,openfpm::grow_policy_identity>(v2);
			}
		}
	}

	// Receiving size
	openfpm::vector<size_t> recv_sz;

	// Receiving buffer for particles ghost get
	openfpm::vector<HeapMemory> recv_mem_gg;

	/*! \brief Call-back to allocate buffer to receive incoming objects (particles)
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
		vector_dist<point,prop,Box,Decomposition,Memory,with_id> * v = static_cast<vector_dist<point,prop,Box,Decomposition,Memory,with_id> *>(ptr);

		v->recv_sz.resize(v->dec.getNNProcessors());
		v->recv_mem_gg.resize(v->dec.getNNProcessors());

		// Get the local processor id
		size_t lc_id = v->dec.ProctoID(i);

		// resize the receive buffer
		v->recv_mem_gg.get(lc_id).resize(msg_i);
		v->recv_sz.get(lc_id) = msg_i;

		return v->recv_mem_gg.get(lc_id).getPointer();
	}

	// Heap memory receiver
	HeapMemory hp_recv;

	// vector v_proc
	openfpm::vector<size_t> v_proc;

	// Receive counter
	size_t recv_cnt;

	/*! \brief Message allocation
	 *
	 * \param message size required to receive from i
	 * \param total message size to receive from all the processors
	 * \param the total number of processor want to communicate with you
	 * \param i processor id
	 * \param ri request id (it is an id that goes from 0 to total_p, and is unique
	 *           every time message_alloc is called)
	 * \param ptr a pointer to the vector_dist structure
	 *
	 * \return the pointer where to store the message
	 *
	 */
	static void * message_alloc_map(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		// cast the pointer
		vector_dist<point,prop,Box,Decomposition,Memory,with_id> * vd = static_cast<vector_dist<point,prop,Box,Decomposition,Memory,with_id> *>(ptr);

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
	vector_dist_iterator<openfpm::vector<point>> getIterator()
	{
		return vector_dist_iterator<openfpm::vector<point>>(v_pos);
	}

	/*! \brief Get the iterator across the position of the ghost particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<point>> getGhostIterator()
	{
		return vector_dist_iterator<openfpm::vector<point>>(v_pos,g_m);
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
	 * \param File output
	 * \param opt NO_GHOST or WITH_GHOST
	 *
	 * \return if the file has been written correctly
	 *
	 */
	inline bool write(std::string out, int opt = NO_GHOST)
	{
		if (hasEnding(out,".csv"))
		{
			// CSVWriter test
			CSVWriter<openfpm::vector<point>, openfpm::vector<prop> > csv_writer;

			std::string output = std::to_string(v_cl.getProcessUnitID()) + std::string("_") + out;

			// Write the CSV
			return csv_writer.write(output,v_pos.get(INTERNAL),v_prp.get(INTERNAL));
		}

		return false;
	}
};


#endif /* VECTOR_HPP_ */
