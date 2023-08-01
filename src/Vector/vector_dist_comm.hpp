/*
 * vector_dist_comm.hpp
 *
 *  Created on: Aug 18, 2016
 *      Author: i-bird
 */

#ifndef SRC_VECTOR_VECTOR_DIST_COMM_HPP_
#define SRC_VECTOR_VECTOR_DIST_COMM_HPP_

#define TEST1

#if defined(CUDA_GPU) && defined(__NVCC__)
#include "Vector/cuda/vector_dist_cuda_funcs.cuh"
#include "util/cuda/kernels.cuh"
#endif

#include "Vector/util/vector_dist_funcs.hpp"
#include "cuda/vector_dist_comm_util_funcs.cuh"
#include "util/cuda/scan_ofp.cuh"

template<typename T>
struct DEBUG
{
	static float ret(T & tmp)
	{
		return 0.0;
	}
};

template<>
struct DEBUG<float &>
{
	static float ret(float & tmp)
	{
		return tmp;
	}
};

/*! \brief compute the communication options from the ghost_get/put options
 *
 *
 */
inline static size_t compute_options(size_t opt)
{
	size_t opt_ = NONE;
	if (opt & NO_CHANGE_ELEMENTS && opt & SKIP_LABELLING)
	{opt_ = RECEIVE_KNOWN | KNOWN_ELEMENT_OR_BYTE;}

	if (opt & RUN_ON_DEVICE)
	{
#if defined(CUDA_GPU) && defined(__NVCC__)
		// Before doing the communication on RUN_ON_DEVICE we have to be sure that the previous kernels complete
		opt_ |= MPI_GPU_DIRECT;
#else
		std::cout << __FILE__ << ":" << __LINE__ << " error: to use the option RUN_ON_DEVICE you must compile with NVCC" << std::endl;
#endif
	}

	return opt_;
}

/*! \brief template selector for asynchronous or not asynchronous
 *
 * \tparam impl implementation
 * \tparam prp properties
 *
 */
template<unsigned int impl, template<typename> class layout_base, unsigned int ... prp>
struct ghost_exchange_comm_impl
{
	template<typename Vcluster_type, typename vector_prop_type,
	         typename vector_pos_type, typename send_vector,
	         typename prc_recv_get_type, typename prc_g_opart_type,
	         typename recv_sz_get_type, typename recv_sz_get_byte_type,
	         typename g_opart_sz_type>
	static inline void sendrecv_prp(Vcluster_type & v_cl,
						 openfpm::vector<send_vector> & g_send_prp,
						 vector_prop_type & v_prp,
						 vector_pos_type & v_pos,
						 prc_g_opart_type & prc_g_opart,
						 prc_recv_get_type & prc_recv_get,
						 recv_sz_get_type & recv_sz_get,
						 recv_sz_get_byte_type & recv_sz_get_byte,
						 g_opart_sz_type & g_opart_sz,
						 size_t g_m,
						 size_t opt)
	{
		// if there are no properties skip
		// SSendRecvP send everything when we do not give properties

		if (sizeof...(prp) != 0)
		{
			size_t opt_ = compute_options(opt);
			if (opt & SKIP_LABELLING)
			{
				if (opt & RUN_ON_DEVICE)
				{
					op_ssend_gg_recv_merge_run_device opm(g_m);
					v_cl.template SSendRecvP_op<op_ssend_gg_recv_merge_run_device,send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,opm,prc_recv_get,recv_sz_get,opt_);
				}
				else
				{
					op_ssend_gg_recv_merge opm(g_m);
					v_cl.template SSendRecvP_op<op_ssend_gg_recv_merge,send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,opm,prc_recv_get,recv_sz_get,opt_);
				}
			}
			else
			{v_cl.template SSendRecvP<send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,prc_recv_get,recv_sz_get,recv_sz_get_byte,opt_);}

			// fill g_opart_sz
			g_opart_sz.resize(prc_g_opart.size());

			for (size_t i = 0 ; i < prc_g_opart.size() ; i++)
				g_opart_sz.get(i) = g_send_prp.get(i).size();
		}
	}

	template<typename Vcluster_type, typename vector_prop_type,
		     typename vector_pos_type, typename send_pos_vector,
		     typename prc_recv_get_type, typename prc_g_opart_type,
		     typename recv_sz_get_type>
	static inline void sendrecv_pos(Vcluster_type & v_cl,
									openfpm::vector<send_pos_vector> & g_pos_send,
									vector_prop_type & v_prp,
									vector_pos_type & v_pos,
									prc_recv_get_type & prc_recv_get,
									recv_sz_get_type & recv_sz_get,
									prc_g_opart_type & prc_g_opart,
									size_t opt)
	{
		size_t opt_ = compute_options(opt);
		if (opt & SKIP_LABELLING)
		{
			v_cl.template SSendRecv<send_pos_vector,decltype(v_pos),layout_base>(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get,opt_);
		}
		else
		{
			prc_recv_get.clear();
			recv_sz_get.clear();
			v_cl.template SSendRecv<send_pos_vector,decltype(v_pos),layout_base>(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get,opt_);
		}
	}

	template<typename Vcluster_type, typename vector_prop_type,
		     typename vector_pos_type, typename send_pos_vector,
		     typename prc_recv_get_type, typename prc_g_opart_type,
		     typename recv_sz_get_type>
	static inline void sendrecv_pos_wait(Vcluster_type & v_cl,
										 openfpm::vector<send_pos_vector> & g_pos_send,
										 vector_prop_type & v_prp,
										 vector_pos_type & v_pos,
										 prc_recv_get_type & prc_recv_get,
										 recv_sz_get_type & recv_sz_get,
										 prc_g_opart_type & prc_g_opart,
										 size_t opt)
	{}

	template<typename Vcluster_type, typename vector_prop_type,
	         typename vector_pos_type, typename send_vector,
	         typename prc_recv_get_type, typename prc_g_opart_type,
	         typename recv_sz_get_type, typename recv_sz_get_byte_type,
	         typename g_opart_sz_type>
	static inline void sendrecv_prp_wait(Vcluster_type & v_cl,
			 	 	 	 	 	 	 	 openfpm::vector<send_vector> & g_send_prp,
			 	 	 	 	 	 	 	 vector_prop_type & v_prp,
			 	 	 	 	 	 	 	 vector_pos_type & v_pos,
			 	 	 	 	 	 	 	 prc_g_opart_type & prc_g_opart,
			 	 	 	 	 	 	 	 prc_recv_get_type & prc_recv_get,
			 	 	 	 	 	 	 	 recv_sz_get_type & recv_sz_get,
			 	 	 	 	 	 	 	 recv_sz_get_byte_type & recv_sz_get_byte,
			 	 	 	 	 	 	 	 g_opart_sz_type & g_opart_sz,
			 	 	 	 	 	 	 	 size_t g_m,
			 	 	 	 	 	 	 	 size_t opt)
	{}
};


template<template<typename> class layout_base, unsigned int ... prp>
struct ghost_exchange_comm_impl<GHOST_ASYNC,layout_base, prp ... >
{
	template<typename Vcluster_type, typename vector_prop_type,
	         typename vector_pos_type, typename send_vector,
	         typename prc_recv_get_type, typename prc_g_opart_type,
	         typename recv_sz_get_type, typename recv_sz_get_byte_type,
	         typename g_opart_sz_type>
	static inline void sendrecv_prp(Vcluster_type & v_cl,
						 openfpm::vector<send_vector> & g_send_prp,
						 vector_prop_type & v_prp,
						 vector_pos_type & v_pos,
						 prc_g_opart_type & prc_g_opart,
						 prc_recv_get_type & prc_recv_get,
						 recv_sz_get_type & recv_sz_get,
						 recv_sz_get_byte_type & recv_sz_get_byte,
						 g_opart_sz_type & g_opart_sz,
						 size_t g_m,
						 size_t opt)
	{
		prc_recv_get.clear();
		recv_sz_get.clear();

		// if there are no properties skip
		// SSendRecvP send everything when we do not give properties

		if (sizeof...(prp) != 0)
		{
			size_t opt_ = compute_options(opt);
			if (opt & SKIP_LABELLING)
			{
				if (opt & RUN_ON_DEVICE)
				{
					op_ssend_gg_recv_merge_run_device opm(g_m);
					v_cl.template SSendRecvP_opAsync<op_ssend_gg_recv_merge_run_device,send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,opm,prc_recv_get,recv_sz_get,opt_);
				}
				else
				{
					op_ssend_gg_recv_merge opm(g_m);
					v_cl.template SSendRecvP_opAsync<op_ssend_gg_recv_merge,send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,opm,prc_recv_get,recv_sz_get,opt_);
				}
			}
			else
			{v_cl.template SSendRecvPAsync<send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,prc_recv_get,recv_sz_get,recv_sz_get_byte,opt_);}
		}

		// fill g_opart_sz
		g_opart_sz.resize(prc_g_opart.size());

		for (size_t i = 0 ; i < prc_g_opart.size() ; i++)
		{g_opart_sz.get(i) = g_send_prp.get(i).size();}
	}

	template<typename Vcluster_type, typename vector_prop_type,
		     typename vector_pos_type, typename send_pos_vector,
		     typename prc_recv_get_type, typename prc_g_opart_type,
		     typename recv_sz_get_type>
	static inline void sendrecv_pos(Vcluster_type & v_cl,
									openfpm::vector<send_pos_vector> & g_pos_send,
									vector_prop_type & v_prp,
									vector_pos_type & v_pos,
									prc_recv_get_type & prc_recv_get,
									recv_sz_get_type & recv_sz_get,
									prc_g_opart_type & prc_g_opart,
									size_t opt)
	{
		prc_recv_get.clear();
		recv_sz_get.clear();

		size_t opt_ = compute_options(opt);
		if (opt & SKIP_LABELLING)
		{
			v_cl.template SSendRecvAsync<send_pos_vector,decltype(v_pos),layout_base>(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get,opt_);
		}
		else
		{
			prc_recv_get.clear();
			recv_sz_get.clear();
			v_cl.template SSendRecvAsync<send_pos_vector,decltype(v_pos),layout_base>(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get,opt_);
		}
	}

	template<typename Vcluster_type, typename vector_prop_type,
		     typename vector_pos_type, typename send_pos_vector,
		     typename prc_recv_get_type, typename prc_g_opart_type,
		     typename recv_sz_get_type>
	static inline void sendrecv_pos_wait(Vcluster_type & v_cl,
										 openfpm::vector<send_pos_vector> & g_pos_send,
										 vector_prop_type & v_prp,
										 vector_pos_type & v_pos,
										 prc_recv_get_type & prc_recv_get,
										 recv_sz_get_type & recv_sz_get,
										 prc_g_opart_type & prc_g_opart,
										 size_t opt)
	{
		size_t opt_ = compute_options(opt);
		if (opt & SKIP_LABELLING)
		{
			v_cl.template SSendRecvWait<send_pos_vector,decltype(v_pos),layout_base>(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get,opt_);
		}
		else
		{
			v_cl.template SSendRecvWait<send_pos_vector,decltype(v_pos),layout_base>(g_pos_send,v_pos,prc_g_opart,prc_recv_get,recv_sz_get,opt_);
		}
	}

	template<typename Vcluster_type, typename vector_prop_type,
	         typename vector_pos_type, typename send_vector,
	         typename prc_recv_get_type, typename prc_g_opart_type,
	         typename recv_sz_get_type, typename recv_sz_get_byte_type,
	         typename g_opart_sz_type>
	static inline void sendrecv_prp_wait(Vcluster_type & v_cl,
			 	 	 	 	 	 	 	 openfpm::vector<send_vector> & g_send_prp,
			 	 	 	 	 	 	 	 vector_prop_type & v_prp,
			 	 	 	 	 	 	 	 vector_pos_type & v_pos,
			 	 	 	 	 	 	 	 prc_g_opart_type & prc_g_opart,
			 	 	 	 	 	 	 	 prc_recv_get_type & prc_recv_get,
			 	 	 	 	 	 	 	 recv_sz_get_type & recv_sz_get,
			 	 	 	 	 	 	 	 recv_sz_get_byte_type & recv_sz_get_byte,
			 	 	 	 	 	 	 	 g_opart_sz_type & g_opart_sz,
			 	 	 	 	 	 	 	 size_t g_m,
			 	 	 	 	 	 	 	 size_t opt)
	{
		// if there are no properties skip
		// SSendRecvP send everything when we do not give properties

		if (sizeof...(prp) != 0)
		{
			size_t opt_ = compute_options(opt);
			if (opt & SKIP_LABELLING)
			{
				if (opt & RUN_ON_DEVICE)
				{
					op_ssend_gg_recv_merge_run_device opm(g_m);
					v_cl.template SSendRecvP_opWait<op_ssend_gg_recv_merge_run_device,send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,opm,prc_recv_get,recv_sz_get,opt_);
				}
				else
				{
					op_ssend_gg_recv_merge opm(g_m);
					v_cl.template SSendRecvP_opWait<op_ssend_gg_recv_merge,send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,opm,prc_recv_get,recv_sz_get,opt_);
				}
			}
			else
			{v_cl.template SSendRecvPWait<send_vector,decltype(v_prp),layout_base,prp...>(g_send_prp,v_prp,prc_g_opart,prc_recv_get,recv_sz_get,recv_sz_get_byte,opt_);}
		}
	}
};


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

template<unsigned int dim,
         typename St,
         typename prop,
         typename Decomposition = CartDecomposition<dim,St>,
         typename Memory = HeapMemory,
         template<typename> class layout_base = memory_traits_lin>
class vector_dist_comm
{
	//! Number of units for each sub-domain
	size_t v_sub_unit_factor = 64;

	//! definition of the send vector for position
	typedef openfpm::vector<Point<dim, St>,Memory,layout_base,openfpm::grow_policy_identity> send_pos_vector;

	//! VCluster
	Vcluster<Memory> & v_cl;

	//! Domain decomposition
	Decomposition dec;

	//! It map the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! For each near processor, outgoing particle id
	//! \warning opart is assumed to be an ordered list
	//! first id particle id
	//! second id shift id
	//! third id is the processor id
	openfpm::vector<aggregate<int,int,int>,
					Memory,
					layout_base > m_opart;

	//! Per processor ordered particles id for ghost_get (see prc_g_opart)
	//! For each processor the internal vector store the id of the
	//! particles that must be communicated to the other processors
	openfpm::vector<openfpm::vector<aggregate<size_t,size_t>>> g_opart;

	//! Same as g_opart but on device, the vector of vector is flatten into a single vector
    openfpm::vector<aggregate<unsigned int,unsigned long int>,
                    CudaMemory,
                    memory_traits_inte> g_opart_device;

	//! Helper buffer for computation (on GPU) of local particles (position)
	openfpm::vector<Point<dim, St>,Memory,layout_base> v_pos_tmp;

	//! Helper buffer for computation (on GPU) of local particles (properties)
	openfpm::vector<prop,Memory,layout_base> v_prp_tmp;

	//! Per processor number of particle g_opart_sz.get(i) = g_opart.get(i).size()
	openfpm::vector<size_t> g_opart_sz;

	//! processor rank list of g_opart
	openfpm::vector<size_t> prc_g_opart;

	//! It store the list of processor that communicate with us (local processor)
	//! from the last ghost get
	openfpm::vector<size_t> prc_recv_get_pos;
	openfpm::vector<size_t> prc_recv_get_prp;

	//! the same as prc_recv_get but for put
	openfpm::vector<size_t> prc_recv_put;

	//! the same as prc_recv_get but for map
	openfpm::vector<size_t> prc_recv_map;

	//! It store the size of the elements added for each processor that communicate with us (local processor)
	//! from the last ghost get
	openfpm::vector<size_t> recv_sz_get_pos;
	openfpm::vector<size_t> recv_sz_get_prp;
	//! Conversion to byte of recv_sz_get
	openfpm::vector<size_t> recv_sz_get_byte;


	//! The same as recv_sz_get but for put
	openfpm::vector<size_t> recv_sz_put;

	//! The same as recv_sz_get but for map
	openfpm::vector<size_t> recv_sz_map;

	//! elements sent for each processors (ghost_get)
	openfpm::vector<size_t> prc_sz_gg;

	//! temporary buffer to processors ids
    openfpm::vector<aggregate<unsigned int>,
                            Memory,
                            layout_base> proc_id_out;

    //! temporary buffer for the scan result
	openfpm::vector<aggregate<unsigned int>,
                             Memory,
                             layout_base> starts;

	//! Processor communication size
	openfpm::vector<aggregate<unsigned int, unsigned int>,Memory,layout_base> prc_offset;


	//! Temporary CudaMemory to do stuff
	CudaMemory mem;

	//! Local ghost marker (across the ghost particles it mark from where we have the)
	//! replicated ghost particles that are local
	size_t lg_m;

	//! Sending buffer
	openfpm::vector_fr<Memory> hsmem;

	//! process the particle with properties
	template<typename prp_object, int ... prp>
	struct proc_with_prp
	{
		//! process the particle
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

	/*! \brief Get the number of particles received from each processor during the last ghost_get
	 *
	 *
	 * \param i processor (list index)
	 * \return the number of particles
	 */
	size_t get_last_ghost_get_received_parts(size_t i)
	{
		// If the last ghost_get did not have properties the information about the number of particles
		// received is in recv_sz_get_ois
		if (recv_sz_get_prp.size() != 0)
		{return recv_sz_get_prp.get(i);}
		else
		{return recv_sz_get_pos.get(i);}
	}

	/*! \brief Get the number of processor involved during the last ghost_get
	 *
	 * \return the number of processor
	 */
	size_t get_last_ghost_get_num_proc()
	{
		if (prc_recv_get_prp.size() != 0)
		{return prc_recv_get_prp.size();}
		else
		{return prc_recv_get_pos.size();}
	}

	/*! \brief Get the number of processor involved during the last ghost_get
	 *
	 * \return the number of processor
	 */
	openfpm::vector<size_t> & get_last_ghost_get_num_proc_vector()
	{
		if (prc_recv_get_prp.size() != 0)
		{return prc_recv_get_prp;}
		else
		{return prc_recv_get_pos;}
	}

	/*! \brief Calculate sending buffer size for each processor
	 *
	 * \param prc_sz_r processor size
	 * \param prc_r processor ids
	 *
	 */
	inline void calc_send_buffers(openfpm::vector<aggregate<unsigned int,unsigned int>,Memory,layout_base> & prc_sz,
								  openfpm::vector<size_t> & prc_sz_r,
								  openfpm::vector<size_t> & prc_r,
								  size_t opt)
	{
		if (opt & RUN_ON_DEVICE)
		{
#ifndef TEST1
			size_t prev_off = 0;
			for (size_t i = 0; i < prc_sz.size() ; i++)
			{
				if (prc_sz.template get<1>(i) != (unsigned int)-1)
				{
					prc_r.add(prc_sz.template get<1>(i));
					prc_sz_r.add(prc_sz.template get<0>(i) - prev_off);
				}
				prev_off = prc_sz.template get<0>(i);
			}
#else

			// Calculate the sending buffer size for each processor, put this information in
			// a contiguous buffer

			for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
			{
				if (prc_sz.template get<0>(i) != 0 && v_cl.rank() != i)
				{
					prc_r.add(i);
					prc_sz_r.add(prc_sz.template get<0>(i));
				}
			}

#endif
		}
		else
		{
			// Calculate the sending buffer size for each processor, put this information in
			// a contiguous buffer

			p_map_req.resize(v_cl.getProcessingUnits());
			for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
			{
				if (prc_sz.template get<0>(i) != 0)
				{
					p_map_req.get(i) = prc_r.size();
					prc_r.add(i);
					prc_sz_r.add(prc_sz.template get<0>(i));
				}
			}
		}
	}

	//! From which decomposition the shift boxes are calculated
	long int shift_box_ndec = -1;

	//! this map is used to check if a combination is already present
	std::unordered_map<size_t, size_t> map_cmb;

	//! The boxes touching the border of the domain are divided in groups (first vector)
	//! each group contain internal ghost coming from sub-domains of the same section
	openfpm::vector_std<openfpm::vector_std<Box<dim, St>>> box_f;

	//! The boxes touching the border of the domain + shift vector linearized from where they come from
	openfpm::vector<Box<dim, St>,Memory,layout_base> box_f_dev;
	openfpm::vector<aggregate<unsigned int>,Memory,layout_base> box_f_sv;

	//! Store the sector for each group (previous vector)
	openfpm::vector_std<comb<dim>> box_cmb;

	//! Id of the local particle to replicate for ghost_get
	openfpm::vector<aggregate<unsigned int,unsigned int>,Memory,layout_base> o_part_loc;

	//! Processor communication size
	openfpm::vector<aggregate<unsigned int, unsigned int>,Memory,layout_base> prc_sz;

	/*! \brief For every internal ghost box we create a structure that order such internal local ghost box in
	 *         shift vectors
	 *
	 */
	void createShiftBox()
	{
		if (shift_box_ndec == (long int)dec.get_ndec())
		{return;}

		struct sh_box
		{
			size_t shift_id;

			unsigned int box_f_sv;
			Box<dim,St> box_f_dev;

			bool operator<(const sh_box & tmp) const
			{
				return shift_id < tmp.shift_id;
			}

		};
		openfpm::vector<sh_box> reord_shift;
		box_f.clear();
		map_cmb.clear();
		box_cmb.clear();

		// Add local particles coming from periodic boundary, the only boxes that count are the one
		// touching the border
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

				reord_shift.add();
				reord_shift.last().shift_id = dec.getLocalIGhostPos(i, j).lin();
				reord_shift.last().box_f_dev = dec.getLocalIGhostBox(i, j);
				reord_shift.last().box_f_sv = dec.convertShift(dec.getLocalIGhostPos(i, j));
			}
		}

		// now we sort box_f by shift_id, the reason is that we have to avoid duplicated particles
		reord_shift.sort();

		box_f_dev.resize(reord_shift.size());
		box_f_sv.resize(reord_shift.size());

		for (size_t i = 0 ; i < reord_shift.size() ; i++)
		{
			box_f_dev.get(i) = reord_shift.get(i).box_f_dev;
			box_f_sv.template get<0>(i) = reord_shift.get(i).box_f_sv;
		}

#ifdef CUDA_GPU

		// move box_f_dev and box_f_sv to device
		box_f_dev.template hostToDevice<0,1>();
		box_f_sv.template hostToDevice<0>();

#endif

		shift_box_ndec = dec.get_ndec();
	}

	/*! \brief Local ghost from labeled particles
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particles properties
	 * \param opt options
	 *
	 */
	void local_ghost_from_opart(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			                    openfpm::vector<prop,Memory,layout_base> & v_prp,
			                    size_t opt)
	{
		// get the shift vectors
		const openfpm::vector<Point<dim, St>,Memory,layout_base> & shifts = dec.getShiftVectors();

		if (!(opt & NO_POSITION))
		{
			if (opt & RUN_ON_DEVICE)
			{
				local_ghost_from_opart_impl<true,dim,St,prop,Memory,layout_base,std::is_same<Memory,CudaMemory>::value>
				::run(o_part_loc,shifts,v_pos,v_prp,opt);
			}
			else
			{
				for (size_t i = 0 ; i < o_part_loc.size() ; i++)
				{
					size_t lin_id = o_part_loc.template get<1>(i);
					size_t key = o_part_loc.template get<0>(i);

					Point<dim, St> p = v_pos.get(key);
					// shift
					p -= shifts.get(lin_id);

					// add this particle shifting its position
					v_pos.add(p);
					v_prp.get(lg_m+i) = v_prp.get(key);
				}
			}
		}
		else
		{
			if (opt & RUN_ON_DEVICE)
			{
				local_ghost_from_opart_impl<false,dim,St,prop,Memory,layout_base,std::is_same<Memory,CudaMemory>::value>
				::run(o_part_loc,shifts,v_pos,v_prp,opt);
			}
			else
			{
				for (size_t i = 0 ; i < o_part_loc.size() ; i++)
				{
					size_t key = o_part_loc.template get<0>(i);

					v_prp.get(lg_m+i) = v_prp.get(key);
				}
			}
		}
	}

	/*! \brief Local ghost from decomposition
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param g_m ghost marker
	 *
	 */
	void local_ghost_from_dec(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			                  openfpm::vector<prop,Memory,layout_base> & v_prp,
			                  size_t g_m,size_t opt)
	{
		o_part_loc.clear();

		// get the shift vectors
		const openfpm::vector<Point<dim,St>,Memory,layout_base> & shifts = dec.getShiftVectors();

		if (opt & RUN_ON_DEVICE)
		{
			local_ghost_from_dec_impl<dim,St,prop,Memory,layout_base,std::is_same<Memory,CudaMemory>::value>
			::run(o_part_loc,shifts,box_f_dev,box_f_sv,v_cl,starts,v_pos,v_prp,g_m,opt);
		}
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
						if (box_f.get(i).get(j).isInsideNP(v_pos.get(key)) == true)
						{
							size_t lin_id = dec.convertShift(box_cmb.get(i));

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
	void add_loc_particles_bc(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			                  openfpm::vector<prop,Memory,layout_base> & v_prp ,
			                  size_t & g_m,
			                  size_t opt)
	{
		// Create the shift boxes
		createShiftBox();

		if (!(opt & SKIP_LABELLING))
			lg_m = v_prp.size();

		if (box_f.size() == 0)
			return;
		else
		{
			if (opt & SKIP_LABELLING)
			{local_ghost_from_opart(v_pos,v_prp,opt);}
			else
			{local_ghost_from_dec(v_pos,v_prp,g_m,opt);}
		}
	}

	/*! \brief This function fill the send buffer for the particle position after the particles has been label with labelParticles
	 *
	 * \param v_pos vector of particle positions
	 * \param g_pos_send Send buffer to fill
	 *
	 */
	void fill_send_ghost_pos_buf(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
								 openfpm::vector<size_t> & prc_sz,
			                     openfpm::vector<send_pos_vector> & g_pos_send,
			                     size_t opt,
			                     bool async)
	{
		// get the shift vectors
		const openfpm::vector<Point<dim,St>,Memory,layout_base> & shifts = dec.getShiftVectors();

		// create a number of send buffers equal to the near processors
		g_pos_send.resize(prc_sz.size());

		size_t old_hsmem_size = 0;

		// if we do async
		if (async == true)
		{
			old_hsmem_size = hsmem.size();
			resize_retained_buffer(hsmem,g_pos_send.size() + hsmem.size());
		}
		else
		{resize_retained_buffer(hsmem,g_pos_send.size());}

		for (size_t i = 0; i < g_pos_send.size(); i++)
		{
			// Buffer must retained and survive the destruction of the
			// vector
			if (hsmem.get(i+old_hsmem_size).ref() == 0)
			{hsmem.get(i+old_hsmem_size).incRef();}

			// Set the memory for retain the send buffer
			g_pos_send.get(i).setMemory(hsmem.get(i+old_hsmem_size));

			// resize the sending vector (No allocation is produced)
			g_pos_send.get(i).resize(prc_sz.get(i));
		}

		if (opt & RUN_ON_DEVICE)
		{
#if defined(CUDA_GPU) && defined(__NVCC__)

			size_t offset = 0;

			// Fill the sending buffers
			for (size_t i = 0 ; i < g_pos_send.size() ; i++)
			{
				auto ite = g_pos_send.get(i).getGPUIterator();

				CUDA_LAUNCH((process_ghost_particles_pos<dim,decltype(g_opart_device.toKernel()),decltype(g_pos_send.get(i).toKernel()),decltype(v_pos.toKernel()),decltype(shifts.toKernel())>),
				ite,
				g_opart_device.toKernel(), g_pos_send.get(i).toKernel(),
				 v_pos.toKernel(),shifts.toKernel(),offset);

				offset += prc_sz.get(i);
			}

#else

			std::cout << __FILE__ << ":" << __LINE__ << " error RUN_ON_DEVICE require that you compile with NVCC, but it seem compiled with a normal compiler" << std::endl;

#endif
		}
		else
		{
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
	}

	/*! \brief This function fill the send buffer for ghost_put
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param v_prp vector of particle properties
	 * \param g_send_prp Send buffer to fill
	 * \param g_m ghost marker
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp>
	void fill_send_ghost_put_prp_buf(openfpm::vector<prop,Memory,layout_base> & v_prp,
									 openfpm::vector<send_vector> & g_send_prp,
									 size_t & g_m,
									 size_t opt)
	{
		// create a number of send buffers equal to the near processors
		// from which we received

		// NOTE in some case the information can be in prc_recv_get_pos

		size_t nproc = get_last_ghost_get_num_proc();

		g_send_prp.resize(nproc);

		resize_retained_buffer(hsmem,g_send_prp.size());

		for (size_t i = 0; i < g_send_prp.size(); i++)
		{
			// Buffer must retained and survive the destruction of the
			// vector
			if (hsmem.get(i).ref() == 0)
				hsmem.get(i).incRef();

			// Set the memory for retain the send buffer
			g_send_prp.get(i).setMemory(hsmem.get(i));

			size_t n_part_recv = get_last_ghost_get_received_parts(i);

			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(n_part_recv);
		}

		size_t accum = g_m;

		if (opt & RUN_ON_DEVICE)
		{
#if defined(CUDA_GPU) && defined(__NVCC__)

				if (sizeof...(prp) != 0)
				{
					// Fill the sending buffers
					for (size_t i = 0 ; i < g_send_prp.size() ; i++)
					{
						size_t n_part_recv = get_last_ghost_get_received_parts(i);

						auto ite = g_send_prp.get(i).getGPUIterator();

						if (ite.nblocks() == 0) {continue;}

						CUDA_LAUNCH((process_ghost_particles_prp_put<decltype(g_send_prp.get(i).toKernel()),decltype(v_prp.toKernel()),prp...>),
						ite,
						g_send_prp.get(i).toKernel(),
						v_prp.toKernel(),accum);

						accum = accum + n_part_recv;
					}
				}

#else

				std::cout << __FILE__ << ":" << __LINE__ << " error RUN_ON_DEVICE require that you compile with NVCC, but it seem compiled with a normal compiler" << std::endl;

#endif
		}
		else
		{
			// Fill the send buffer
			for (size_t i = 0; i < g_send_prp.size(); i++)
			{
				size_t j2 = 0;
				size_t n_part_recv = get_last_ghost_get_received_parts(i);

				for (size_t j = accum; j < accum + n_part_recv; j++)
				{
					// source object type
					typedef encapc<1, prop, typename openfpm::vector<prop,Memory,layout_base>::layout_type> encap_src;
					// destination object type
					typedef encapc<1, prp_object, typename openfpm::vector<prp_object,Memory,layout_base>::layout_type> encap_dst;

					// Copy only the selected properties
					object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(j), g_send_prp.get(i).get(j2));

					j2++;
				}

				accum = accum + n_part_recv;
			}
		}
	}

	/*! \brief resize the retained buffer by nbf
	 *
	 *
	 */
	void resize_retained_buffer(openfpm::vector_fr<Memory> & rt_buf, size_t nbf)
	{
		// Release all the buffer that are going to be deleted
		for (size_t i = nbf ; i < rt_buf.size() ; i++)
		{
			rt_buf.get(i).decRef();
		}

		hsmem.resize(nbf);
	}

	/*! \brief Set the buffer for each property
	 *
	 *
	 */
	template<typename send_vector, typename v_mpl>
	struct set_mem_retained_buffers_inte
	{
		openfpm::vector<send_vector> & g_send_prp;

		size_t i;

		openfpm::vector_fr<Memory> & hsmem;

		size_t j;

		set_mem_retained_buffers_inte(openfpm::vector<send_vector> & g_send_prp, size_t i ,
				                      openfpm::vector_fr<Memory> & hsmem, size_t j)
		:g_send_prp(g_send_prp),i(i),hsmem(hsmem),j(j)
		{}

		//! It call the setMemory function for each property
		template<typename T>
		inline void operator()(T& t)
		{
			g_send_prp.get(i).template setMemory<T::value>(hsmem.get(j));

			j++;
		}
	};

	template<bool inte_or_lin,typename send_vector, typename v_mpl>
	struct set_mem_retained_buffers
	{
		static inline size_t set_mem_retained_buffers_(openfpm::vector<send_vector> & g_send_prp,
				     	 	 	 	 	 	 openfpm::vector<size_t> & prc_sz,
											 size_t i,
											 openfpm::vector_fr<Memory> & hsmem,
											 size_t j)
		{
			// Set the memory for retain the send buffer
			g_send_prp.get(i).setMemory(hsmem.get(j));

			// resize the sending vector (No allocation is produced)
			g_send_prp.get(i).resize(prc_sz.get(i));

			return j+1;
		}
	};

	template<typename send_vector, typename v_mpl>
	struct set_mem_retained_buffers<true,send_vector,v_mpl>
	{
		static inline size_t set_mem_retained_buffers_(openfpm::vector<send_vector> & g_send_prp,
											 openfpm::vector<size_t> & prc_sz,
				 	 	 	 	 	 	 	 size_t i,
				 	 	 	 	 	 	 	 openfpm::vector_fr<Memory> & hsmem,
				 	 	 	 	 	 	 	 size_t j)
		{
			set_mem_retained_buffers_inte<send_vector,v_mpl> smrbi(g_send_prp,i,hsmem,j);

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,boost::mpl::size<v_mpl>::type::value>>(smrbi);

			// if we do not send properties do not reallocate
			if (boost::mpl::size<v_mpl>::type::value != 0)
			{
				// resize the sending vector (No allocation is produced)
				g_send_prp.get(i).resize(prc_sz.get(i));
			}

			return smrbi.j;
		}
	};

	/*! \brief This function fill the send buffer for properties after the particles has been label with labelParticles
	 *
	 * \tparam send_vector type used to send data
	 * \tparam prp_object object containing only the properties to send
	 * \tparam prp set of properties to send
	 *
	 * \param v_prp vector of particle properties
	 * \param g_send_prp Send buffer to fill
	 *
	 */
	template<typename send_vector, typename prp_object, int ... prp>
	void fill_send_ghost_prp_buf(openfpm::vector<prop,Memory,layout_base> & v_prp,
								 openfpm::vector<size_t> & prc_sz,
			                     openfpm::vector<send_vector> & g_send_prp,
			                     size_t opt)
	{
		size_t factor = 1;

		typedef typename to_boost_vmpl<prp...>::type v_mpl;

		if (is_layout_inte<layout_base<prop>>::value == true) {factor *= sizeof...(prp);}

		// create a number of send buffers equal to the near processors
		g_send_prp.resize(prc_sz.size());

		resize_retained_buffer(hsmem,g_send_prp.size()*factor);

		for (size_t i = 0; i < hsmem.size(); i++)
		{
			// Buffer must retained and survive the destruction of the
			// vector
			if (hsmem.get(i).ref() == 0)
			{hsmem.get(i).incRef();}
		}

		size_t j = 0;
		for (size_t i = 0; i < g_send_prp.size(); i++)
		{
			j = set_mem_retained_buffers<is_layout_inte<layout_base<prop>>::value,send_vector,v_mpl>::set_mem_retained_buffers_(g_send_prp,prc_sz,i,hsmem,j);
		}

		if (opt & RUN_ON_DEVICE)
		{
#if defined(CUDA_GPU) && defined(__NVCC__)

			size_t offset = 0;

			if (sizeof...(prp) != 0)
			{
				// Fill the sending buffers
				for (size_t i = 0 ; i < g_send_prp.size() ; i++)
				{
					auto ite = g_send_prp.get(i).getGPUIterator();

					CUDA_LAUNCH((process_ghost_particles_prp<decltype(g_opart_device.toKernel()),decltype(g_send_prp.get(i).toKernel()),decltype(v_prp.toKernel()),prp...>),
					ite,
					g_opart_device.toKernel(), g_send_prp.get(i).toKernel(),
					 v_prp.toKernel(),offset);

					offset += prc_sz.get(i);
				}
			}

#else

			std::cout << __FILE__ << ":" << __LINE__ << " error RUN_ON_DEVICE require that you compile with NVCC, but it seem compiled with a normal compiler" << std::endl;

#endif
		}
		else
		{
			// if no properties must be sent skip this step
			if (sizeof...(prp) == 0)	{return;}

			// Fill the send buffer
			for (size_t i = 0; i < g_opart.size(); i++)
			{
				for (size_t j = 0; j < g_opart.get(i).size(); j++)
				{
					// source object type
					typedef decltype(v_prp.get(g_opart.get(i).template get<0>(j))) encap_src;
					// destination object type
					typedef decltype(g_send_prp.get(i).get(j)) encap_dst;

					// Copy only the selected properties
					object_si_d<encap_src, encap_dst, OBJ_ENCAP, prp...>(v_prp.get(g_opart.get(i).template get<0>(j)), g_send_prp.get(i).get(j));
				}
			}
		}
	}

	/*! \brief allocate and fill the send buffer for the map function
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particles properties
	 * \param prc_sz_r For each processor in the list the size of the message to send
	 * \param m_pos sending buffer for position
	 * \param m_prp sending buffer for properties
	 * \param offset from where start the list of the particles that migrate in o_part
	 *        This parameter is used only in case of RUN_ON_DEVICE option
	 *
	 */
	void fill_send_map_buf(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			               openfpm::vector<prop,Memory,layout_base> & v_prp,
			               openfpm::vector<size_t> & prc_sz_r,
			               openfpm::vector<size_t> & prc_r,
			               openfpm::vector<openfpm::vector<Point<dim,St>,Memory,layout_base,openfpm::grow_policy_identity>> & m_pos,
			               openfpm::vector<openfpm::vector<prop,Memory,layout_base,openfpm::grow_policy_identity>> & m_prp,
			               openfpm::vector<aggregate<unsigned int, unsigned int>,Memory,layout_base> & prc_sz,
			               size_t opt)
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

		if (opt & RUN_ON_DEVICE)
		{
			if (v_cl.size() == 1)
			{return;}

#if defined(CUDA_GPU) && defined(__NVCC__)

			// The first part of m_opart and prc_sz contain the local particles

			int rank = v_cl.rank();

			v_pos_tmp.resize(prc_sz.template get<0>(rank));
			v_prp_tmp.resize(prc_sz.template get<0>(rank));

			auto ite = v_pos_tmp.getGPUIterator();

			starts.template deviceToHost<0>();
			size_t offset = starts.template get<0>(rank);

			// no work to do
			if (ite.wthr.x != 0)
			{
				// fill v_pos_tmp and v_prp_tmp with local particles
				CUDA_LAUNCH((process_map_particles<decltype(m_opart.toKernel()),decltype(v_pos_tmp.toKernel()),decltype(v_prp_tmp.toKernel()),
					                                           decltype(v_pos.toKernel()),decltype(v_prp.toKernel())>),
				ite,
				m_opart.toKernel(),v_pos_tmp.toKernel(), v_prp_tmp.toKernel(),
					            v_pos.toKernel(),v_prp.toKernel(),offset);
			}

			// Fill the sending buffers
			for (size_t i = 0 ; i < m_pos.size() ; i++)
			{
				size_t offset = starts.template get<0>(prc_r.template get<0>(i));

				auto ite = m_pos.get(i).getGPUIterator();

				// no work to do
				if (ite.wthr.x != 0)
				{

					CUDA_LAUNCH((process_map_particles<decltype(m_opart.toKernel()),decltype(m_pos.get(i).toKernel()),decltype(m_prp.get(i).toKernel()),
						                                           decltype(v_pos.toKernel()),decltype(v_prp.toKernel())>),
					ite,
					m_opart.toKernel(),m_pos.get(i).toKernel(), m_prp.get(i).toKernel(),
						            v_pos.toKernel(),v_prp.toKernel(),offset);

				}
			}

			// old local particles with the actual local particles
			v_pos_tmp.swap(v_pos);
			v_prp_tmp.swap(v_prp);

#else

			std::cout << __FILE__ << ":" << __LINE__ << " error RUN_ON_DEVICE require that you compile with NVCC, but it seem compiled with a normal compiler" << std::endl;

#endif
		}
		else
		{
			// end vector point
			long int id_end = v_pos.size();

			// end opart point
			long int end = m_opart.size()-1;

			// Run through all the particles and fill the sending buffer
			for (size_t i = 0; i < m_opart.size(); i++)
			{
				process_map_particle<proc_without_prp>(i,end,id_end,m_opart,p_map_req,m_pos,m_prp,v_pos,v_prp,cnt);
			}

			v_pos.resize(v_pos.size() - m_opart.size());
			v_prp.resize(v_prp.size() - m_opart.size());
		}
	}


	/*! \brief allocate and fill the send buffer for the map function
	 *
	 * \tparam prp_object object type to send
	 * \tparam prp properties to send
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param prc_sz_r number of particles to send for each processor
	 * \param m_pos sending buffer for position
	 * \param m_prp sending buffer for properties
	 *
	 */
	template<typename prp_object,int ... prp>
	void fill_send_map_buf_list(openfpm::vector<Point<dim, St>> & v_pos,
			                    openfpm::vector<prop,Memory,layout_base> & v_prp,
								openfpm::vector<size_t> & prc_sz_r,
								openfpm::vector<openfpm::vector<Point<dim,St>>> & m_pos,
								openfpm::vector<openfpm::vector<prp_object>> & m_prp)
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
			process_map_particle<proc_with_prp<prp_object,prp...>>(i,end,id_end,m_opart,p_map_req,m_pos,m_prp,v_pos,v_prp,cnt);
		}

		v_pos.resize(v_pos.size() - m_opart.size());
		v_prp.resize(v_prp.size() - m_opart.size());
	}

	/*! \brief Label particles for mappings
	 *
	 * \param v_pos vector of particle positions
	 * \param lbl_p Particle labeled
	 * \param prc_sz For each processor the number of particles to send
	 * \param opt options
	 *
	 */
	template<typename obp> void labelParticleProcessor(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			                                           openfpm::vector<aggregate<int,int,int>,
			                                                           Memory,
			                                                           layout_base> & lbl_p,
			                                           openfpm::vector<aggregate<unsigned int,unsigned int>,Memory,layout_base> & prc_sz,
			                                           size_t opt)
	{
		if (opt == RUN_ON_DEVICE)
		{
#ifdef __NVCC__

			// Map directly on gpu

			lbl_p.resize(v_pos.size());

			// labelling kernel

			prc_sz.template fill<0>(0);

			auto ite = v_pos.getGPUIterator();
			if (ite.wthr.x == 0)
			{
				starts.resize(v_cl.size());
				starts.template fill<0>(0);
				return;
			}

			// we have one process we can skip ...
			if (v_cl.size() == 1)
			{
				// ... but we have to apply the boundary conditions

				periodicity_int<dim> bc;

				for (size_t i = 0 ; i < dim ; i++)	{bc.bc[i] = dec.periodicity(i);}

				CUDA_LAUNCH((apply_bc_each_part<dim,St,decltype(v_pos.toKernel())>),ite,dec.getDomain(),bc,v_pos.toKernel());

				return;
			}

			// label particle processor
			CUDA_LAUNCH((process_id_proc_each_part<dim,St,decltype(dec.toKernel()),decltype(v_pos.toKernel()),decltype(lbl_p.toKernel()),decltype(prc_sz.toKernel())>),
			ite,
			dec.toKernel(),v_pos.toKernel(),lbl_p.toKernel(),prc_sz.toKernel(),v_cl.rank());

			starts.resize(v_cl.size());
			openfpm::scan((unsigned int *)prc_sz.template getDeviceBuffer<0>(), prc_sz.size(), (unsigned int *)starts.template getDeviceBuffer<0>() , v_cl.getgpuContext());

			// move prc_sz to host
			prc_sz.template deviceToHost<0>();

			ite = lbl_p.getGPUIterator();

			// we order lbl_p
			CUDA_LAUNCH((reorder_lbl<decltype(lbl_p.toKernel()),decltype(starts.toKernel())>),ite,lbl_p.toKernel(),starts.toKernel());


#else

			std::cout << __FILE__ << ":" << __LINE__ << " error, it seems you tried to call map with RUN_ON_DEVICE option, this requires to compile the program with NVCC" << std::endl;

#endif
		}
		else
		{
			// reset lbl_p
			lbl_p.clear();
			prc_sz_gg.clear();
			o_part_loc.clear();
			g_opart.clear();
			prc_g_opart.clear();

			// resize the label buffer
			prc_sz.template fill<0>(0);

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
				{p_id = dec.processorID(v_pos.get(key));}
				else
				{p_id = obp::out(key, v_cl.getProcessUnitID());}

				// Particle to move
				if (p_id != v_cl.getProcessUnitID())
				{
					if ((long int) p_id != -1)
					{
						prc_sz.template get<0>(p_id)++;
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
	}

	/*! \brief Label the particles
	 *
	 * It count the number of particle to send to each processors and save its ids
	 *
	 * \see nn_prcs::getShiftvectors()
	 *
	 * \param v_pos vector of particle positions
	 * \param v_prp vector of particle properties
	 * \param prc for each particle it label the processor id (the owner of the particle, or where it should go the particle)
	 * \param g_m ghost marker
	 * \param opt ghost_get options
	 *
	 */
	void labelParticlesGhost(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			                 openfpm::vector<prop,Memory,layout_base> & v_prp,
			                 openfpm::vector<size_t> & prc,
			                 openfpm::vector<size_t> & prc_sz,
			                 openfpm::vector<aggregate<unsigned int,unsigned int>,Memory,layout_base> & prc_offset,
			                 size_t & g_m,
			                 size_t opt)
	{
		// Buffer that contain for each processor the id of the particle to send
		prc_sz.clear();
		g_opart.clear();
		g_opart.resize(dec.getNNProcessors());
		prc_g_opart.clear();

		if (opt & RUN_ON_DEVICE)
		{
			labelParticlesGhost_impl<dim,St,prop,Memory,layout_base,
			                         Decomposition,std::is_same<Memory,CudaMemory>::value>
			::run(mem,dec,g_opart_device,proc_id_out,starts,v_cl,v_pos,v_prp,prc,prc_sz,prc_offset,g_m,opt);
		}
		else
		{
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
					prc_sz.add(g_opart.get(i).size());
					g_opart_f.add();
					g_opart.get(i).swap(g_opart_f.last());
					prc.add(dec.IDtoProc(i));
				}
			}

			g_opart.swap(g_opart_f);
		}
#ifdef EXTREA_TRACE_PRE_COMM
		Extrae_user_function (0);
#endif
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
		vector_dist_comm<dim, St, prop, Decomposition, Memory, layout_base> * vd = static_cast<vector_dist_comm<dim, St, prop, Decomposition, Memory, layout_base> *>(ptr);

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
	vector_dist_comm(const vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> & v)
	:v_cl(create_vcluster<Memory>()),dec(create_vcluster()),lg_m(0)
	{
		this->operator=(v);
	}


	/*! \brief Constructor
	 *
	 * \param dec Domain decompositon
	 *
	 */
	vector_dist_comm(const Decomposition & dec)
	:v_cl(create_vcluster<Memory>()),dec(dec),lg_m(0)
	{

	}

	/*! \brief Constructor
	 *
	 * \param dec Domain decompositon
	 *
	 */
	vector_dist_comm(Decomposition && dec)
	:v_cl(create_vcluster<Memory>()),dec(dec),lg_m(0)
	{

	}

	/*! \brief Constructor
	 *
	 */
	vector_dist_comm()
	:v_cl(create_vcluster<Memory>()),dec(create_vcluster()),lg_m(0)
	{
	}

	/*! \brief Destructor
	 *
	 * Release the retained buffer
	 *
	 */
	~vector_dist_comm()
	{
		for (size_t i = 0 ; i < hsmem.size() ; i++)
		{
			if (hsmem.get(i).ref() == 1)
				hsmem.get(i).decRef();
			else
				std::cout << __FILE__ << ":" << __LINE__ << " internal error memory is in an invalid state " << std::endl;
		}

	}

	/*! \brief Get the number of minimum sub-domain per processor
	 *
	 * \return minimum number
	 *
	 */
	size_t getDecompositionGranularity()
	{
		return v_sub_unit_factor;
	}

	/*! \brief Set the minimum number of sub-domain per processor
	 *
	 * \param n_sub
	 *
	 */
	void setDecompositionGranularity(size_t n_sub)
	{
		this->v_sub_unit_factor = n_sub;
	}

	/*! \brief Initialize the decomposition
	 *
	 * \param box domain
	 * \param bc boundary conditions
	 * \param g ghost extension
	 * \param opt additional options
	 *
	 */
	void init_decomposition(Box<dim,St> & box,
							const size_t (& bc)[dim],
							const Ghost<dim,St> & g,
							size_t opt,
							const grid_sm<dim,void> & gdist)
	{
		size_t div[dim];

		if (opt & BIND_DEC_TO_GHOST)
		{
			// padding
			size_t pad = 0;

			// CellDecomposer
			CellDecomposer_sm<dim,St,shift<dim,St>> cd_sm;

			// Calculate the divisions for the symmetric Cell-lists
			cl_param_calculateSym<dim,St>(box,cd_sm,g,pad);

			for (size_t i = 0 ; i < dim ; i++)
			{div[i] = cd_sm.getDiv()[i] - 2*pad;}

			// Create the sub-domains
			dec.setParameters(div, box, bc, g, gdist);
		}
		else
		{
			dec.setGoodParameters(box, bc, g, getDecompositionGranularity(), gdist);
		}
		dec.decompose();
	}

	/*! \brief Initialize the decomposition
	 *
	 * \param box domain
	 * \param bc boundary conditions
	 * \param g ghost extension
	 * \param opt additional options
	 *
	 */
	void init_decomposition_gr_cell(Box<dim,St> & box,
							const size_t (& bc)[dim],
							const Ghost<dim,St> & g,
							size_t opt,
							const grid_sm<dim,void> & gdist)
	{
		size_t div[dim];

		for (size_t i = 0 ; i < dim ; i++)
		{div[i] = gdist.size(i);}

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
	template<unsigned int impl, int ... prp> inline void ghost_get_(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
												 openfpm::vector<prop,Memory,layout_base> & v_prp,
												 size_t & g_m,
												 size_t opt = WITH_POSITION)
	{
#ifdef PROFILE_SCOREP
		SCOREP_USER_REGION("ghost_get",SCOREP_USER_REGION_TYPE_FUNCTION)
#endif

		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object,Memory,layout_base,openfpm::grow_policy_identity> send_vector;

		if (!(opt & NO_POSITION))
		{v_pos.resize(g_m);}

		// reset the ghost part

		if (!(opt & SKIP_LABELLING))
		{v_prp.resize(g_m);}

		// Label all the particles
		if ((opt & SKIP_LABELLING) == false)
		{labelParticlesGhost(v_pos,v_prp,prc_g_opart,prc_sz_gg,prc_offset,g_m,opt);}

		{
			// Send and receive ghost particle information
			openfpm::vector<send_vector> g_send_prp;

			fill_send_ghost_prp_buf<send_vector, prp_object, prp...>(v_prp,prc_sz_gg,g_send_prp,opt);

	#if defined(CUDA_GPU) && defined(__NVCC__)
			cudaDeviceSynchronize();
	#endif

			// if there are no properties skip
			// SSendRecvP send everything when we do not give properties

			ghost_exchange_comm_impl<impl,layout_base,prp ...>::template
			sendrecv_prp(v_cl,g_send_prp,v_prp,v_pos,prc_g_opart,
					 prc_recv_get_prp,recv_sz_get_prp,recv_sz_get_byte,g_opart_sz,g_m,opt);
		}

		if (!(opt & NO_POSITION))
		{
			// Sending buffer for the ghost particles position
			openfpm::vector<send_pos_vector> g_pos_send;

			fill_send_ghost_pos_buf(v_pos,prc_sz_gg,g_pos_send,opt,impl == GHOST_ASYNC);

#if defined(CUDA_GPU) && defined(__NVCC__)
			cudaDeviceSynchronize();
#endif

			ghost_exchange_comm_impl<impl,layout_base,prp ...>::template
			sendrecv_pos(v_cl,g_pos_send,v_prp,v_pos,prc_recv_get_pos,recv_sz_get_pos,prc_g_opart,opt);

            // fill g_opart_sz
            g_opart_sz.resize(prc_g_opart.size());

			for (size_t i = 0 ; i < prc_g_opart.size() ; i++)
				g_opart_sz.get(i) = g_pos_send.get(i).size();
		}

        // Important to ensure that the number of particles in v_prp must be equal to v_pos
        // Note that if we do not give properties sizeof...(prp) == 0 in general at this point
        // v_prp.size() != v_pos.size()
        if (!(opt & SKIP_LABELLING))
        {
                v_prp.resize(v_pos.size());
        }

		add_loc_particles_bc(v_pos,v_prp,g_m,opt);
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
	template<int ... prp> inline void ghost_wait_(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 openfpm::vector<prop,Memory,layout_base> & v_prp,
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 size_t & g_m,
			 	 	 	 	 	 	 	 	 	 	 	 	 	 	 size_t opt = WITH_POSITION)
	{
		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object,Memory,layout_base,openfpm::grow_policy_identity> send_vector;

		// Send and receive ghost particle information
		openfpm::vector<send_vector> g_send_prp;
		openfpm::vector<send_pos_vector> g_pos_send;

		ghost_exchange_comm_impl<GHOST_ASYNC,layout_base,prp ...>::template
		sendrecv_prp_wait(v_cl,g_send_prp,v_prp,v_pos,prc_g_opart,
					 prc_recv_get_prp,recv_sz_get_prp,recv_sz_get_byte,g_opart_sz,g_m,opt);


		ghost_exchange_comm_impl<GHOST_ASYNC,layout_base,prp ...>::template
		sendrecv_pos_wait(v_cl,g_pos_send,v_prp,v_pos,prc_recv_get_pos,recv_sz_get_pos,prc_g_opart,opt);
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
	 * \param opt options
	 *
	 */
	template<unsigned int ... prp> void map_list_(openfpm::vector<Point<dim, St>> & v_pos, openfpm::vector<prop> & v_prp, size_t & g_m, size_t opt)
	{
		if (opt & RUN_ON_DEVICE)
		{
			std::cout << "Error: " << __FILE__ << ":" << __LINE__ << " map_list is unsupported on device (coming soon)" << std::endl;
			return;
		}

		typedef KillParticle obp;

		// Processor communication size
		openfpm::vector<aggregate<unsigned int,unsigned int>,Memory,layout_base> prc_sz(v_cl.getProcessingUnits());

		// map completely reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// m_opart, Contain the processor id of each particle (basically where they have to go)
		labelParticleProcessor<obp>(v_pos,m_opart, prc_sz,opt);

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		p_map_req.resize(v_cl.getProcessingUnits());
		openfpm::vector<size_t> prc_sz_r;
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (prc_sz.template get<0>(i) != 0)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(prc_sz.template get<0>(i));
			}
		}

		if (opt & MAP_LOCAL)
		{
			// if the map is local we indicate that we receive only from the neighborhood processors

			prc_recv_map.clear();
			for (size_t i = 0 ; i < dec.getNNProcessors() ; i++)
			{prc_recv_map.add(dec.IDtoProc(i));}
		}

		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		//! position vector
		openfpm::vector<openfpm::vector<Point<dim, St>>> m_pos;
		//! properties vector
		openfpm::vector<openfpm::vector<prp_object>> m_prp;

		fill_send_map_buf_list<prp_object,prp...>(v_pos,v_prp,prc_sz_r, m_pos, m_prp);

		v_cl.SSendRecv(m_pos,v_pos,prc_r,prc_recv_map,recv_sz_map,opt);
		v_cl.template SSendRecvP<openfpm::vector<prp_object>,decltype(v_prp),layout_base,prp...>(m_prp,v_prp,prc_r,prc_recv_map,recv_sz_map,opt);

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
	template<typename obp = KillParticle>
	void map_(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
			  openfpm::vector<prop,Memory,layout_base> & v_prp, size_t & g_m,
			  size_t opt)
	{
#ifdef PROFILE_SCOREP
		SCOREP_USER_REGION("map",SCOREP_USER_REGION_TYPE_FUNCTION)
#endif

		prc_sz.resize(v_cl.getProcessingUnits());

		// map completely reset the ghost part
		v_pos.resize(g_m);
		v_prp.resize(g_m);

		// Contain the processor id of each particle (basically where they have to go)
		labelParticleProcessor<obp>(v_pos,m_opart, prc_sz,opt);

		openfpm::vector<size_t> prc_sz_r;
		openfpm::vector<size_t> prc_r;

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		calc_send_buffers(prc_sz,prc_sz_r,prc_r,opt);

		//! position vector
		openfpm::vector<openfpm::vector<Point<dim, St>,Memory,layout_base,openfpm::grow_policy_identity>> m_pos;
		//! properties vector
		openfpm::vector<openfpm::vector<prop,Memory,layout_base,openfpm::grow_policy_identity>> m_prp;

		fill_send_map_buf(v_pos,v_prp, prc_sz_r,prc_r, m_pos, m_prp,prc_sz,opt);

		size_t opt_ = 0;
		if (opt & RUN_ON_DEVICE)
		{
#if defined(CUDA_GPU) && defined(__NVCC__)
			// Before doing the communication on RUN_ON_DEVICE we have to be sure that the previous kernels complete
			cudaDeviceSynchronize();
			opt_ |= MPI_GPU_DIRECT;
#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: to use the option RUN_ON_DEVICE you must compile with NVCC" << std::endl;
#endif
		}

		v_cl.template SSendRecv<openfpm::vector<Point<dim, St>,Memory,layout_base,openfpm::grow_policy_identity>,
					   openfpm::vector<Point<dim, St>,Memory,layout_base>,
					   layout_base>
					   (m_pos,v_pos,prc_r,prc_recv_map,recv_sz_map,opt_);

		v_cl.template SSendRecv<openfpm::vector<prop,Memory,layout_base,openfpm::grow_policy_identity>,
					   openfpm::vector<prop,Memory,layout_base>,
					   layout_base>
					   (m_prp,v_prp,prc_r,prc_recv_map,recv_sz_map,opt_);

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
	vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> & operator=(const vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> & vc)
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
	vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> & operator=(vector_dist_comm<dim,St,prop,Decomposition,Memory,layout_base> && vc)
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
	 * \param opt options
	 *
	 */
	template<template<typename,typename> class op, int ... prp>
	void ghost_put_(openfpm::vector<Point<dim, St>,Memory,layout_base> & v_pos,
					openfpm::vector<prop,Memory,layout_base> & v_prp,
					size_t & g_m,
					size_t opt)
	{
		// Sending property object
		typedef object<typename object_creator<typename prop::type, prp...>::type> prp_object;

		// send vector for each processor
		typedef openfpm::vector<prp_object,Memory,layout_base> send_vector;

		openfpm::vector<send_vector> g_send_prp;
		fill_send_ghost_put_prp_buf<send_vector, prp_object, prp...>(v_prp,g_send_prp,g_m,opt);

		if (opt & RUN_ON_DEVICE)
		{
#if defined(CUDA_GPU) && defined(__NVCC__)
			// Before doing the communication on RUN_ON_DEVICE we have to be sure that the previous kernels complete
			cudaDeviceSynchronize();
#else
			std::cout << __FILE__ << ":" << __LINE__ << " error: to use the option RUN_ON_DEVICE you must compile with NVCC" << std::endl;
#endif
		}

		// Send and receive ghost particle information
		if (opt & NO_CHANGE_ELEMENTS)
		{
			size_t opt_ = compute_options(opt);

			if (opt & RUN_ON_DEVICE)
			{
				op_ssend_recv_merge_gpu<op,decltype(g_opart_device),decltype(prc_offset)> opm(g_opart_device,prc_offset);
				v_cl.template SSendRecvP_op<op_ssend_recv_merge_gpu<op,decltype(g_opart_device),decltype(prc_offset)>,
				                            send_vector,
											decltype(v_prp),
											layout_base,
											prp...>(g_send_prp,v_prp,prc_recv_get_prp,opm,prc_g_opart,g_opart_sz,opt_);
			}
			else
			{
				op_ssend_recv_merge<op,decltype(g_opart)> opm(g_opart);
				v_cl.template SSendRecvP_op<op_ssend_recv_merge<op,decltype(g_opart)>,
				                            send_vector,
											decltype(v_prp),
											layout_base,
											prp...>(g_send_prp,v_prp,prc_recv_get_prp,opm,prc_g_opart,g_opart_sz,opt_);
			}
		}
		else
		{
			size_t opt_ = compute_options(opt);

			if (opt & RUN_ON_DEVICE)
			{
				op_ssend_recv_merge_gpu<op,decltype(g_opart_device),decltype(prc_offset)> opm(g_opart_device,prc_offset);
				v_cl.template SSendRecvP_op<op_ssend_recv_merge_gpu<op,decltype(g_opart_device),decltype(prc_offset)>,
				                            send_vector,
											decltype(v_prp),
											layout_base,
											prp...>(g_send_prp,v_prp,get_last_ghost_get_num_proc_vector(),opm,prc_recv_put,recv_sz_put,opt_);
			}
			else
			{
				op_ssend_recv_merge<op,decltype(g_opart)> opm(g_opart);
				v_cl.template SSendRecvP_op<op_ssend_recv_merge<op,decltype(g_opart)>,
				                            send_vector,
											decltype(v_prp),
											layout_base,
											prp...>(g_send_prp,v_prp,get_last_ghost_get_num_proc_vector(),opm,prc_recv_put,recv_sz_put,opt_);
			}
		}

		// process also the local replicated particles

		if (lg_m < v_prp.size() && v_prp.size() - lg_m != o_part_loc.size())
		{
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Local ghost particles = " << v_prp.size() - lg_m << " != " << o_part_loc.size() << std::endl;
			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Check that you did a ghost_get before a ghost_put" << std::endl;
		}


		if (opt & RUN_ON_DEVICE)
		{
			v_prp.template merge_prp_v_device<op,prop,Memory,
											  openfpm::grow_policy_double,
											  layout_base,
											  decltype(o_part_loc),prp ...>(v_prp,lg_m,o_part_loc);
		}
		else
		{
			v_prp.template merge_prp_v<op,prop,Memory,
			                           openfpm::grow_policy_double,
									   layout_base,
									   decltype(o_part_loc),prp ...>(v_prp,lg_m,o_part_loc);
		}
	}
};


#endif /* SRC_VECTOR_VECTOR_DIST_COMM_HPP_ */
