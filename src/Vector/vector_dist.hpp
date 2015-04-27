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
#include "common.hpp"

#define NO_ID false
#define ID true

#define GET	1
#define PUT 2

/*! \brief Distributed vector
 *
 */

template<typename point, typename prop, typename Box, typename Decomposition , typename Memory=HeapMemory, bool with_id=false>
class vector_dist
{
private:

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
	 * \param number of elements
	 *
	 */
	vector_dist(size_t np, Box box)
	:dec(Decomposition(*global_v_cluster)),v_cl(*global_v_cluster)
	{
		typedef ::Box<point::dims,typename point::coord_type> b;

		// Allocate unassigned particles vectors
		v_pos = v_cl.template allocate<openfpm::vector<point>>(1);
		v_prp = v_cl.template allocate<openfpm::vector<prop>>(1);

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
		for (int i = 0 ; i < point::dims ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/point::dims));}

		// Create the sub-domains
		dec.setParameters(div,box);

		// Get the bounding box containing the processor domain +- one sub-domain spacing
		::Box<point::dims,typename point::coord_type> & bbound = dec.getProcessorBounds();

		// the smallest sub-division of the domain on each dimension
		typename point::coord_type smallest_doms[point::dims];

		// convert spacing divisions
		size_t n_g[point::dims];

		for (size_t i = 0 ; i < point::dims ; i++)
		{
			n_g[i] = box.template getBase<b::p2>(i) / smallest_doms[i];
		}

		point p;
		p.zero();

		// Initialize the geo cell list
		geo_cell.Initialize(box,n_g,p,8);
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
	template<unsigned int id> inline auto getProp(vect_dist_key_dx vec_key) -> decltype(v_prp.get(vec_key.v_c).template get<id>(vec_key.key))
	{
		return v_prp.get(vec_key.v_c).template get<id>(vec_key.key);
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

		// Create the sz and prc buffer

		openfpm::vector<size_t> prc_sz_r;
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0 ; i < v_cl.getProcessingUnits() ; i++)
		{
			if (p_map.get(i) == 1)
			{
				prc_r.add(i);
				prc_sz_r.add(prc_sz.get(i));
			}
		}

		// Allocate all the buffers

		openfpm::vector<pos_prop> pb(prc_r.size());

		for (size_t i = 0 ;  i < prc_r.size() ; i++)
		{
			// Create the size required to store the particles position and properties to communicate
			size_t s1 = openfpm::vector<point>::calculateMem(prc_sz_r.get(i),0);
			size_t s2 = openfpm::vector<prop>::calculateMem(prc_sz_r.get(i),0);

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

			lbl = (lbl > v_cl.getProcessUnitID())?lbl-1:lbl;

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

		for (size_t i = 0 ; i < v_cl.getProcessingUnits() ; i++)
		{
			prc_sz_r.get(i) = prc_sz_r.get(i)*(sizeof(prop) + sizeof(point));
		}

		// Send and receive the particles

		recv_cnt = 0;
		v_cl.sendrecvMultipleMessages(prc_sz_r.size(),&p_map.get(0), &prc_sz_r.get(0), &prc_r.get(0) , &ptr.get(0) , vector_dist::message_alloc, this ,NEED_ALL_SIZE);

		// overwrite the outcoming particle with the incoming particle and resize the vectors

		size_t o_p_id = 0;

		for (size_t i = 0 ; i < v_proc.size() ; i++)
		{
			// Get the number of elements

			size_t n_ele = v_proc.get(i) / (sizeof(point) + sizeof(prop));

			PtrMemory * ptr1 = new PtrMemory(hp_recv.getPointer(),n_ele * sizeof(point));
			PtrMemory * ptr2 = new PtrMemory((unsigned char *)hp_recv.getPointer() + n_ele * sizeof(point),n_ele * sizeof(prop));

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
		}

		// remove the hole (out-going particles) in the vector

		v_pos.get(0).remove(opart,o_p_id);
		v_prp.get(0).remove(opart,o_p_id);
	}

	// ghost particles sending buffer
	openfpm::vector<HeapMemory> ghost_send_hp;

	// Each entry contain the size of the ghost sending buffer
	std::unordered_map<size_t,size_t> ghost_prc_sz;

	// ghost particle labels
	openfpm::vector<size_t> ghost_lbl_p;

	/*! \brief It synchronize getting the ghost particles
	 *
	 * \prp Properties to get
	 * \opt options
	 * 		NO_RELABEL: If the particles does not move avoid to relabel
	 *
	 */
	template<unsigned int N> void ghost_get(const size_t prp[N], size_t opt)
	{
		// outgoing particles-id
		openfpm::vector<size_t> opart;

		ghost_prc_sz.clear();

		// Label the internal (assigned) particles
		auto it = v_pos.get(0).getIterator();

		// Label all the particles with the processor id, where they should go
		while (it.isNext())
		{
			auto key = it.get();

			size_t p_id = dec.ghost_processorID(v_pos.get(0).get(key));

			ghost_lbl_p.get(key) = p_id;

			// It has to communicate
			if (p_id != v_cl.getProcessUnitID())
			{
				// add particle to communicate
				ghost_prc_sz[p_id]++;

				opart.add(key);
			}

			++it;
		}

		// Create memory allocator for the send buffers
		size_t i = 0;
		ghost_send_hp.resize(ghost_prc_sz.size());

		for ( auto it = ghost_prc_sz.begin(); it != ghost_prc_sz.end(); ++it )
		{
			// we are sending only some properties, so calculate the size of the sending buffer
			size_t element_size = ele_size<N,typename prop::type>(prp);

			ghost_send_hp.get(i).resize(it->second * element_size);

			i++;
		}

		//

		// ca

		// send and receive the properties of the particles

		// add the received particles to the vector
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
	 * \param ptr a pointer to the vector_dist structure
	 *
	 * \return the pointer where to store the message
	 *
	 */
	static void * message_alloc(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, void * ptr)
	{
		// cast the pointer
		vector_dist<point,prop,Box,Decomposition,Memory,with_id> * vd = static_cast<vector_dist<point,prop,Box,Decomposition,Memory,with_id> *>(ptr);

		// Resize the memory and
		vd->hp_recv.resize(total_msg);
		vd->v_proc.resize(total_p);

		// Return the receive pointer
		void * recv_ptr = (unsigned char *)vd->hp_recv.getPointer() + vd->recv_cnt;

		// increment the receive pointer
		vd->recv_cnt += msg_i;

		// Save the processor message size
		vd->v_proc.get(i) = msg_i;

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

	/*! \brief Get the iterator across the properties of the particles
	 *
	 * \return an iterator
	 *
	 */
	vector_dist_iterator<openfpm::vector<point>> getPropIterator()
	{
		return vector_dist_iterator<openfpm::vector<prop>>(v_prp);
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
};


#endif /* VECTOR_HPP_ */
