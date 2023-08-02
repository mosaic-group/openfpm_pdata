/*
 * grid_dist_id_comm.hpp
 *
 *  Created on: Nov 13, 2016
 *      Author: yaroslav
 */

#ifndef SRC_GRID_GRID_DIST_ID_COMM_HPP_
#define SRC_GRID_GRID_DIST_ID_COMM_HPP_

#include "Vector/vector_dist_ofb.hpp"
#include "Grid/copy_grid_fast.hpp"
#include "grid_dist_util.hpp"
#include "util/common_pdata.hpp"
#include "lib/pdata.hpp"
#include "Grid/grid_common.hpp"


/*! \brief Unpack selector
 *
 *
 */
template<bool result,typename T, typename device_grid, typename Memory>
struct grid_unpack_selector_with_prp
{
	/*! \brief Error i do not know how to unpack
	 *
	 * \param recv_buf buffer with data
	 * \param sub2 where to unpack (extension)
	 * \param gd grid where to unpack
	 * \param ps unpack status
	 *
	 */
	template<template<typename,typename> class op, typename sub_it_type, int ... prp> static void call_unpack(ExtPreAlloc<Memory> & recv_buf, sub_it_type & sub2, device_grid & gd, Unpack_stat & ps)
	{
		std::cerr << __FILE__ << ":" << __LINE__ << " Error: complex properties on grids are not supported yet" << std::endl;
	}
};

/*! \brief Unpack selector
 *
 *
 */
template<typename T, typename device_grid, typename Memory>
struct grid_unpack_selector_with_prp<true,T,device_grid,Memory>
{

	/*! \brief Unpack
	 *
	 * \param recv_buf buffer with data
	 * \param sub2 where to unpack (extension)
	 * \param gd grid where to unpack
	 * \param ps unpack status
	 *
	 */
	template<template<typename,typename> class op, typename sub_it_type, unsigned int ... prp>
	static void call_unpack(ExtPreAlloc<Memory> & recv_buf,
							sub_it_type & sub2,
							device_grid & gd,
							Unpack_stat & ps)
	{
		gd.template unpack_with_op<op,Memory,prp ...>(recv_buf,sub2,ps);
	}
};

/*! \brief Unpack selector
 *
 * Stub version
 *
 */
template<typename device_grid, typename Memory, typename T>
struct grid_call_serialize_variadic {};

/*! \brief Unpack selector
 *
 * Selector when there is not max_prop
 *
 */
template<typename device_grid, typename Memory , int ... prp>
struct grid_call_serialize_variadic<device_grid, Memory, index_tuple<prp...>>
{

	/*! \brief Unpack
	 *
	 * \param recv_buf buffer with data
	 * \param sub2 where to unpack (extension)
	 * \param dg grid where to unpack
	 * \param ps unpack status
	 *
	 */
	template<template<typename,typename> class op, typename sub_it_type, typename T>
	inline static void call_unpack(ExtPreAlloc<Memory> & recv_buf,
									sub_it_type & sub2,
									device_grid & dg,
									Unpack_stat & ps)
	{
		const bool result = has_pack_gen<typename T::type>::value == false;

		grid_unpack_selector_with_prp<result,T,device_grid,Memory>::template call_unpack<op,sub_it_type,prp...>(recv_buf,sub2,dg,ps);
	}
};

/*! \brief Unpack selector
 *
 * Selector when there is max_prop
 *
 */
template<template<typename,typename> class op, typename T, typename device_grid, typename Memory>
struct grid_unpack_with_prp
{

	/*! \brief Unpack
	 *
	 * \param recv_buf buffer with data
	 * \param sub2 where to unpack (extension)
	 * \param dg grid where to unpack
	 * \param ps unpack status
	 *
	 */
	template<typename sub_it_type, unsigned int ... prp> static void unpacking(ExtPreAlloc<Memory> & recv_buf, sub_it_type & sub2, device_grid & dg, Unpack_stat & ps)
	{
		typedef index_tuple<prp...> ind_prop_to_pack;
		grid_call_serialize_variadic<device_grid,Memory,ind_prop_to_pack>::template call_unpack<op,sub_it_type,T>(recv_buf, sub2, dg, ps);
	}
};

/*! \brief This class is an helper for the communication of grid_dist_id
 *
 * \tparam dim Dimensionality of the grid
 * \tparam St Type of space where the grid is living
 * \tparam T object the grid is storing
 * \tparam Decomposition Class that decompose the grid for example CartDecomposition
 * \tparam Memory Is the allocator
 * \tparam device_grid of base structure is going to store the data
 *
 * \see grid_dist_id
 *
 */

template<unsigned int dim, typename St, typename T, typename Decomposition = CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_id_comm
{
	//! VCluster
	Vcluster<Memory> & v_cl;

	//! Maps the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! Stores the list of processors that communicate with us (local processor)
	openfpm::vector<size_t> prc_recv_map;

	//! Stores the size of the elements added for each processor that communicate with us (local processor)
	openfpm::vector<size_t> recv_sz_map;

	//! List of processor to send to
	openfpm::vector<size_t> send_prc_queue;

	//! Pointer to the memory to send
	openfpm::vector<void *> send_pointer;

	//! size to send
	openfpm::vector<size_t> send_size;

	//! receiving buffers in case of dynamic
	openfpm::vector_fr<BMemory<Memory>> recv_buffers;

	struct rp_id
	{
		int p_id;
		int size;
		int i;

		bool operator<(const rp_id & tmp) const
		{
			return p_id < tmp.p_id;
		}
	};

	//! receiving processors
	openfpm::vector<rp_id> recv_proc;

	//! For each near processor, outgoing intersection grid
	//! \warning m_oGrid is assumed to be an ordered list
	//! first id is grid
	//! second id is the processor id
	openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> m_oGrid;
	openfpm::vector<int> m_oGrid_c;

	//! Memory for the ghost sending buffer
	Memory g_send_prp_mem;

	//! Memory for the ghost receiving buffer
	Memory g_recv_prp_mem;

	//! send pointers
	openfpm::vector<void *> pointers;
	openfpm::vector<void *> pointers2;

	//! header unpacker info
	openfpm::vector_gpu<aggregate<void *,void *,int>> pointers_h;
	int n_headers_slot = 1;
	openfpm::vector_gpu<aggregate<size_t,size_t,unsigned int>> headers;

	//! Receiving option
	size_t opt;

	/*! \brief Sync the local ghost part
	 *
	 * \tparam prp... properties to sync
	 *
	 * \param loc_ig_box local internel ghost boxes
	 * \param loc_eg_box local external ghost boxes
	 * \param gdb_ext information about the local grids
	 * \param loc_grid local grids
	 * \param g_id_to_external_ghost_box from global index to external ghost box
	 *
	 */
	template<int... prp> void ghost_get_local(const openfpm::vector<i_lbox_grid<dim>> & loc_ig_box,
											  const openfpm::vector<e_lbox_grid<dim>> & loc_eg_box,
											  const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
											  openfpm::vector<device_grid> & loc_grid,
											  std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box,
											  const grid_sm<dim,void> & ginfo,
											  bool use_bx_def,
											  size_t opt)
	{
		rem_copy_opt opt_ = rem_copy_opt::NONE_OPT;
		if (opt & SKIP_LABELLING)
		{opt_ = rem_copy_opt::KEEP_GEOMETRY;}

		if (opt_ != rem_copy_opt::KEEP_GEOMETRY)
		{
			for (size_t i = 0 ; i < loc_grid.size() ; i++)
			{loc_grid.get(i).copyRemoveReset();}
		}

		grid_key_dx<dim> cnt[1];
		cnt[0].zero();

		//! For all the sub-domains
		for (size_t i = 0 ; i < loc_ig_box.size() ; i++)
		{
			//! For all the internal ghost boxes of each sub-domain
			for (size_t j = 0 ; j < loc_ig_box.get(i).bid.size() ; j++)
			{
				size_t sub_id_src_gdb_ext = loc_ig_box.get(i).bid.get(j).sub_gdb_ext;

				// sub domain connected with external box
				size_t sub_id_dst = loc_ig_box.get(i).bid.get(j).sub;

				// local internal ghost box connected
				for (size_t v = 0 ; v < loc_ig_box.get(i).bid.get(j).k.size() ; v++)
				{
					size_t k = loc_ig_box.get(i).bid.get(j).k.get(v);

					Box<dim,long int> bx_dst = loc_eg_box.get(sub_id_dst).bid.get(k).ebox;

					// convert into local
					size_t sub_id_dst_gdb_ext = loc_eg_box.get(sub_id_dst).bid.get(k).sub_gdb_ext;
					bx_dst -= gdb_ext.get(sub_id_dst_gdb_ext).origin;

					// create 2 sub grid iterator

					if (bx_dst.isValid() == false)
					{continue;}

					Box<dim,long int>  bx_src = flip_box(loc_eg_box.get(sub_id_dst).bid.get(k).ebox,loc_eg_box.get(sub_id_dst).bid.get(k).cmb,ginfo);
					bx_src -= gdb_ext.get(sub_id_src_gdb_ext).origin;

	#ifdef SE_CLASS1

					if (use_bx_def == false)
					{
						if (loc_eg_box.get(sub_id_dst).bid.get(k).sub != i)
						{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination are not correctly linked" << "\n";}
					}

					if (bx_src.getVolumeKey() != bx_dst.getVolumeKey())
					{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination does not match in size" << "\n";}

	#endif

					auto & gd = loc_grid.get(sub_id_dst_gdb_ext);

					gd.remove(bx_dst);
					gd.copy_to(loc_grid.get(sub_id_src_gdb_ext),bx_src,bx_dst);
				}
			}
		}

		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{
			loc_grid.get(i).template removeCopyToFinalize<prp ...>(v_cl.getgpuContext(), rem_copy_opt::PHASE1 | opt_);
		}

		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{
			loc_grid.get(i).template removeCopyToFinalize<prp ...>(v_cl.getgpuContext(), rem_copy_opt::PHASE2 | opt_);
		}

		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{
			loc_grid.get(i).template removeCopyToFinalize<prp ...>(v_cl.getgpuContext(), rem_copy_opt::PHASE3 | opt_);
		}
	}

	/*! \brief Sync the local ghost part
	 *
	 * \tparam prp... properties to sync
	 *
	 * \param loc_ig_box local internel ghost boxes
	 * \param loc_eg_box local external ghost boxes
	 * \param gdb_ext information about the local grids
	 * \param loc_grid local grids
	 * \param g_id_to_external_ghost_box global-if to external ghost box
	 *
	 */
	template<template<typename,typename> class op, int... prp> void ghost_put_local(const openfpm::vector<i_lbox_grid<dim>> & loc_ig_box,
											  const openfpm::vector<e_lbox_grid<dim>> & loc_eg_box,
											  const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
											  openfpm::vector<device_grid> & loc_grid,
											  openfpm::vector<std::unordered_map<size_t,size_t>> & g_id_to_external_ghost_box)
	{
		//! For all the sub-domains
		for (size_t i = 0 ; i < loc_eg_box.size() ; i++)
		{
			//! For all the external ghost boxes of each sub-domain
			for (size_t j = 0 ; j < loc_eg_box.get(i).bid.size() ; j++)
			{
				if (loc_eg_box.get(i).bid.get(j).initialized == false)
					continue;

				Box<dim,long int> bx_src = loc_eg_box.get(i).bid.get(j).ebox;
				// convert into local
				bx_src -= gdb_ext.get(i).origin;

				// sub domain connected with external box
				size_t sub_id_dst = loc_eg_box.get(i).bid.get(j).sub;

				// local external ghost box connected
				size_t k = loc_eg_box.get(i).bid.get(j).k;

				Box<dim,long int> bx_dst = loc_ig_box.get(sub_id_dst).bid.get(k).box;

				// convert into local
				bx_dst -= gdb_ext.get(sub_id_dst).origin;

				// create 2 sub grid iterator

				if (bx_dst.isValid() == false)
				{continue;}

#ifdef SE_CLASS1

				if (loc_ig_box.get(sub_id_dst).bid.get(k).sub != i)
					std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination are not correctly linked" << "\n";

				if (bx_src.getVolume() != bx_dst.getVolume())
				{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination does not match in size" << "\n";}

#endif

				auto & gd2 = loc_grid.get(sub_id_dst);
				gd2.template copy_to_op<op,prp...>(loc_grid.get(i),bx_src,bx_dst);

			}
		}
	}

	/* Send or queue the the information
	 *
	 * This function send or queue the information to the other processor. In case the
	 * device grid is a compressed format like in multi-resolution the communication is
	 * queued because the other side does not know the size of the communication. If is
	 * not compressed the other side know the size so a direct send is done
	 *
	 */
	void send_or_queue(size_t prc, char * pointer, char * pointer2)
	{
		if (device_grid::isCompressed() == false)
		{v_cl.send(prc,0,pointer,(char *)pointer2 - (char *)pointer);}
		else
		{
			send_prc_queue.add(prc);
			send_pointer.add(pointer);
			send_size.add(pointer2-pointer);
		}
	}

	static void * receive_dynamic(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, size_t tag, void * ptr)
	{
		grid_dist_id_comm * gd = static_cast<grid_dist_id_comm *>(ptr);

		gd->recv_buffers.add();

		gd->recv_buffers.last().resize(msg_i);
		gd->recv_proc.add();
		gd->recv_proc.last().p_id = i;
		gd->recv_proc.last().size = msg_i;
		gd->recv_proc.last().i = gd->recv_proc.size()-1;

		if (gd->opt & RUN_ON_DEVICE)
		{
			return gd->recv_buffers.last().getDevicePointer();
		}

		return gd->recv_buffers.last().getPointer();
	}

	/* Send or queue the the information
	 *
	 * This function send or queue the information to the other processor. In case the
	 * device grid is a compressed format like in multi-resolution the communication is
	 * queued because the other side does not know the size of the communication. If is
	 * not compressed the other side know the size so a direct send is done
	 *
	 */
	template <typename prp_object>
	void queue_recv_data_get(const openfpm::vector<ep_box_grid<dim>> & eg_box,
			    		 std::vector<size_t> & prp_recv,
						 ExtPreAlloc<Memory> & prRecv_prp)
	{
#ifdef __NVCC__
		cudaDeviceSynchronize();
#endif

		if (device_grid::isCompressed() == false)
		{
			//! Receive the information from each processors
			for ( size_t i = 0 ; i < eg_box.size() ; i++ )
			{
				prp_recv.push_back(eg_box.get(i).recv_pnt * sizeof(prp_object) + sizeof(size_t)*eg_box.get(i).n_r_box);
			}

			size_t tot_recv = ExtPreAlloc<Memory>::calculateMem(prp_recv);

			//! Resize the receiving buffer
			g_recv_prp_mem.resize(tot_recv);

			// queue the receives
			for ( size_t i = 0 ; i < eg_box.size() ; i++ )
			{
				prRecv_prp.allocate(prp_recv[i]);
				v_cl.recv(eg_box.get(i).prc,0,prRecv_prp.getPointer(),prp_recv[i]);
			}
		}
		else
		{
			// It is not possible to calculate the total information so we have to receive

			if (send_prc_queue.size() == 0)
			{
                v_cl.sendrecvMultipleMessagesNBX(send_prc_queue.size(),NULL,
                                                                        NULL,NULL,
                                                                        receive_dynamic,this);
			}
			else
			{
				v_cl.sendrecvMultipleMessagesNBX(send_prc_queue.size(),&send_size.get(0),
											 &send_prc_queue.get(0),&send_pointer.get(0),
											 receive_dynamic,this);
			}

			// Reorder what we received

			recv_proc.sort();

			openfpm::vector_fr<BMemory<Memory>> tmp;
			tmp.resize(recv_proc.size());

			for (int i = 0 ; i < recv_proc.size() ; i++)
			{
				tmp.get(i).swap(recv_buffers.get(recv_proc.get(i).i));
			}

			recv_buffers.swap(tmp);
		}
	}

	/* Send or queue the the information
	 *
	 * This function send or queue the information to the other processor. In case the
	 * device grid is a compressed format like in multi-resolution the communication is
	 * queued because the other side does not know the size of the communication. If is
	 * not compressed the other side know the size so a direct send is done
	 *
	 */
	template <typename prp_object>
	void queue_recv_data_put(const openfpm::vector<ip_box_grid<dim>> & ig_box,
			    		 std::vector<size_t> & prp_recv,
						 ExtPreAlloc<Memory> & prRecv_prp)
	{
		if (device_grid::isCompressed() == false)
		{
			// Receive the information from each processors
			for ( size_t i = 0 ; i < ig_box.size() ; i++ )
			{
				prp_recv.push_back(0);

				// for each external ghost box
				for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
				{
					// External ghost box
					Box<dim,size_t> g_ig_box = ig_box.get(i).bid.get(j).box;
					prp_recv[prp_recv.size()-1] += g_ig_box.getVolumeKey() * sizeof(prp_object) + sizeof(size_t);
				}
			}

			size_t tot_recv = ExtPreAlloc<Memory>::calculateMem(prp_recv);

			//! Resize the receiving buffer
			g_recv_prp_mem.resize(tot_recv);

			prRecv_prp.incRef();

			// queue the receives
			for ( size_t i = 0 ; i < ig_box.size() ; i++ )
			{
				prRecv_prp.allocate(prp_recv[i]);
				v_cl.recv(ig_box.get(i).prc,0,prRecv_prp.getPointer(),prp_recv[i]);
			}

			prRecv_prp.decRef();
		}
		else
		{
			// It is not possible to calculate the total information so we have to receive

			if (send_prc_queue.size() == 0)
			{
				v_cl.sendrecvMultipleMessagesNBX(send_prc_queue.size(),NULL,
											 NULL,NULL,
											 receive_dynamic,this);
			}
			else
			{
				v_cl.sendrecvMultipleMessagesNBX(send_prc_queue.size(),&send_size.get(0),
											 &send_prc_queue.get(0),&send_pointer.get(0),
											 receive_dynamic,this);
			}
		}
	}

	template<typename mem,unsigned ... prp>
	void unpack_data_to_ext_ghost(ExtPreAlloc<mem> & emem,
									openfpm::vector<device_grid> & loc_grid,
									size_t i,
									const openfpm::vector<ep_box_grid<dim>> & eg_box,
									const std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box,
									const openfpm::vector<e_box_multi<dim>> & eb_gid_list,
									Unpack_stat & ps,
									size_t opt)
	{
		// Unpack the ghost box global-id

		size_t g_id;
		// we move from device to host the gid
		if (opt & RUN_ON_DEVICE)
		{emem.deviceToHost(ps.getOffset(),ps.getOffset()+sizeof(size_t));}
		Unpacker<size_t,mem>::unpack(emem,g_id,ps);

		size_t l_id = 0;
		// convert the global id into local id
		auto key = g_id_to_external_ghost_box.find(g_id);

		if (key != g_id_to_external_ghost_box.end()) // FOUND
		{l_id = key->second;}
		else
		{
			// NOT FOUND

			// It must be always found, if not it mean that the processor has no-idea of
			// what is stored and conseguently do not know how to unpack, print a critical error
			// and return

			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Critical, cannot unpack object, because received data cannot be interpreted\n";

			return;
		}


		// we unpack into the last eb_gid_list that is always big enought to
		// unpack the information

		size_t le_id = eb_gid_list.get(l_id).full_match;
		size_t ei =	eb_gid_list.get(l_id).e_id;

		// Get the external ghost box associated with the packed information
		Box<dim,long int> box = eg_box.get(ei).bid.get(le_id).l_e_box;
		size_t sub_id = eg_box.get(ei).bid.get(le_id).sub;

		// sub-grid where to unpack
		auto sub2 = loc_grid.get(sub_id).getIterator(box.getKP1(),box.getKP2(),false);

		rem_copy_opt opt_ = rem_copy_opt::NONE_OPT;
		if (opt & SKIP_LABELLING)
		{opt_ = rem_copy_opt::KEEP_GEOMETRY;}

		// Unpack
		loc_grid.get(sub_id).remove(box);
		Unpacker<device_grid,mem>::template unpack<decltype(sub2),decltype(v_cl.getgpuContext()),prp...>(emem,sub2,loc_grid.get(sub_id),ps,v_cl.getgpuContext(),opt_);

		// Copy the information on the other grid
		for (long int j = 0 ; j < (long int)eb_gid_list.get(l_id).eb_list.size() ; j++)
		{
			size_t nle_id = eb_gid_list.get(l_id).eb_list.get(j);
			if (nle_id != le_id)
			{
//				size_t nle_id = eb_gid_list.get(l_id).eb_list.get(j);
				size_t n_sub_id = eg_box.get(ei).bid.get(nle_id).sub;

				Box<dim,long int> box = eg_box.get(ei).bid.get(nle_id).l_e_box;
				Box<dim,long int> rbox = eg_box.get(ei).bid.get(nle_id).lr_e_box;

				loc_grid.get(n_sub_id).remove(box);
				loc_grid.get(n_sub_id).copy_to(loc_grid.get(sub_id),rbox,box);
			}
		}
	}

	template<typename mem, typename header_type,unsigned ... prp>
	void unpack_data_to_ext_ghost_with_header(ExtPreAlloc<mem> & emem,
									openfpm::vector<device_grid> & loc_grid,
									header_type & headers,
									size_t i,
									const openfpm::vector<ep_box_grid<dim>> & eg_box,
									const std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box,
									const openfpm::vector<e_box_multi<dim>> & eb_gid_list,
									Unpack_stat & ps,
									size_t opt)
	{
		// Unpack the ghost box global-id

		size_t g_id;
		// we move from device to host the gid
		g_id = headers.template get<0>(i);
		ps.addOffset(sizeof(size_t));

		size_t l_id = 0;
		// convert the global id into local id
		auto key = g_id_to_external_ghost_box.find(g_id);

		if (key != g_id_to_external_ghost_box.end()) // FOUND
		{l_id = key->second;}
		else
		{
			// NOT FOUND

			// It must be always found, if not it mean that the processor has no-idea of
			// what is stored and conseguently do not know how to unpack, print a critical error
			// and return

			std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Critical, cannot unpack object, because received data cannot be interpreted\n";

			return;
		}


		// we unpack into the last eb_gid_list that is always big enought to
		// unpack the information

		size_t le_id = eb_gid_list.get(l_id).full_match;
		size_t ei =	eb_gid_list.get(l_id).e_id;

		// Get the external ghost box associated with the packed information
		Box<dim,long int> box = eg_box.get(ei).bid.get(le_id).l_e_box;
		size_t sub_id = eg_box.get(ei).bid.get(le_id).sub;

		// sub-grid where to unpack
		auto sub2 = loc_grid.get(sub_id).getIterator(box.getKP1(),box.getKP2(),false);

		rem_copy_opt opt_ = rem_copy_opt::NONE_OPT;
		if (opt & SKIP_LABELLING)
		{opt_ = rem_copy_opt::KEEP_GEOMETRY;}

		// Unpack
		loc_grid.get(sub_id).remove(box);
		Unpacker<device_grid,mem>::template unpack_with_header<decltype(sub2),decltype(headers),decltype(v_cl.getgpuContext()),prp...>
																				(emem,
																				sub2,
																				loc_grid.get(sub_id),
																				headers,
																				i,
																				ps,
																				v_cl.getgpuContext(),
																				opt_);

		// Copy the information on the other grid
		for (long int j = 0 ; j < (long int)eb_gid_list.get(l_id).eb_list.size() ; j++)
		{
			size_t nle_id = eb_gid_list.get(l_id).eb_list.get(j);
			if (nle_id != le_id)
			{
//				size_t nle_id = eb_gid_list.get(l_id).eb_list.get(j);
				size_t n_sub_id = eg_box.get(ei).bid.get(nle_id).sub;

				Box<dim,long int> box = eg_box.get(ei).bid.get(nle_id).l_e_box;
				Box<dim,long int> rbox = eg_box.get(ei).bid.get(nle_id).lr_e_box;

				loc_grid.get(n_sub_id).remove(box);
				loc_grid.get(n_sub_id).copy_to(loc_grid.get(sub_id),rbox,box);
			}
		}
	}

	template<unsigned int ... prp>
	void fill_headers(size_t opt)
	{
		if ((opt & KEEP_PROPERTIES) == 0 && device_grid::is_unpack_header_supported())
		{
			headers.resize(n_headers_slot * recv_buffers.size());

			Memory result;
			result.allocate(sizeof(int));

			pointers_h.resize(recv_buffers.size());

			for ( size_t i = 0 ; i < recv_buffers.size() ; i++ )
			{
				pointers_h.template get<0>(i) = recv_buffers.get(i).getDevicePointer();
				pointers_h.template get<1>(i) = (unsigned char *)recv_buffers.get(i).getDevicePointer() + recv_buffers.get(i).size();
			}

			pointers_h.template hostToDevice<0,1>();

			while(1)
			{
				for ( size_t i = 0 ; i < recv_buffers.size() ; i++ )
				{pointers_h.template get<2>(i) = 0;}
				pointers_h.template hostToDevice<2>();
				*(int *)result.getPointer() = 0;
				result.hostToDevice();

				device_grid::template unpack_headers<decltype(pointers_h),decltype(headers),decltype(result),prp ...>(pointers_h,headers,result,n_headers_slot);
				result.deviceToHost();

				if (*(int *)result.getPointer() == 0) {break;}

				n_headers_slot *= 2;
				headers.resize(n_headers_slot * recv_buffers.size());

			}

			headers.template deviceToHost<0,1,2>();
		}
	}

	template<unsigned ... prp>
	void merge_received_data_get(openfpm::vector<device_grid> & loc_grid,
							const openfpm::vector<ep_box_grid<dim>> & eg_box,
							const std::vector<size_t> & prp_recv,
							ExtPreAlloc<Memory> & prRecv_prp,
							const std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box,
							const openfpm::vector<e_box_multi<dim>> & eb_gid_list,
							size_t opt)
	{
		if (device_grid::isCompressed() == false)
		{
			// wait to receive communication
			v_cl.execute();

			Unpack_stat ps;

			// Unpack the object
			for ( size_t i = 0 ; i < eg_box.size() ; i++ )
			{
				size_t mark_here = ps.getOffset();

				// for each external ghost box
				while (ps.getOffset() - mark_here < prp_recv[i])
				{
					// Unpack the ghost box global-id


					unpack_data_to_ext_ghost<Memory,prp ...>(prRecv_prp,loc_grid,i,
																eg_box,g_id_to_external_ghost_box,eb_gid_list,
																ps,opt);
				}
			}
		}
		else
		{
			fill_headers<prp ...>(opt);

			if (headers.size() != 0)
			{
				// Unpack the object
				for ( size_t i = 0 ; i < recv_buffers.size() ; i++ )
				{
					Unpack_stat ps;
					size_t mark_here = ps.getOffset();

					ExtPreAlloc<BMemory<Memory>> mem(recv_buffers.get(i).size(),recv_buffers.get(i));

					int j = 0;

					// for each external ghost box
					while (ps.getOffset() - mark_here < recv_buffers.get(i).size())
					{
						// Unpack the ghost box global-id

						unpack_data_to_ext_ghost_with_header<BMemory<Memory>,decltype(headers),prp ...>(mem,loc_grid,headers,i*n_headers_slot+j,
																	eg_box,g_id_to_external_ghost_box,eb_gid_list,
																	ps,opt);

						j++;
					}
				}
			}
			else
			{
				// Unpack the object
				for ( size_t i = 0 ; i < recv_buffers.size() ; i++ )
				{
					Unpack_stat ps;
					size_t mark_here = ps.getOffset();

					ExtPreAlloc<BMemory<Memory>> mem(recv_buffers.get(i).size(),recv_buffers.get(i));

					// for each external ghost box
					while (ps.getOffset() - mark_here < recv_buffers.get(i).size())
					{
						// Unpack the ghost box global-id

						unpack_data_to_ext_ghost<BMemory<Memory>,prp ...>(mem,loc_grid,i,
																	eg_box,g_id_to_external_ghost_box,eb_gid_list,
																	ps,opt);
					}
				}
			}
		}
	}


	template<template<typename,typename> class op, unsigned ... prp>
	void merge_received_data_put(Decomposition & dec, openfpm::vector<device_grid> & loc_grid,
							const openfpm::vector<ip_box_grid<dim>> & ig_box,
							const std::vector<size_t> & prp_recv,
							ExtPreAlloc<Memory> & prRecv_prp,
							const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
							const openfpm::vector<std::unordered_map<size_t,size_t>> & g_id_to_internal_ghost_box)
	{
		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		if (device_grid::isCompressed() == false)
		{
			v_cl.execute();

			Unpack_stat ps;

			// Unpack the object
			for ( size_t i = 0 ; i < ig_box.size() ; i++ )
			{
				// for each external ghost box
				for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
				{
					// Unpack the ghost box global-id

					size_t g_id;
					Unpacker<size_t,HeapMemory>::unpack(prRecv_prp,g_id,ps);

					size_t l_id = 0;
					// convert the global id into local id
					auto key = g_id_to_internal_ghost_box.get(i).find(g_id);
					if (key != g_id_to_internal_ghost_box.get(i).end()) // FOUND
					{l_id = key->second;}
					else
					{
						// NOT FOUND

						// It must be always found, if not it mean that the processor has no-idea of
						// what is stored and conseguently do not know how to unpack, print a critical error
						// and return

						std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Critical, cannot unpack object, because received data cannot be interpreted\n";

						return;
					}

					// Get the internal ghost box associated with the packed information
					Box<dim,size_t> box = ig_box.get(i).bid.get(l_id).box;
					size_t sub_id = ig_box.get(i).bid.get(l_id).sub;
					box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

					// sub-grid where to unpack
					auto sub2 = loc_grid.get(sub_id).getIterator(box.getKP1(),box.getKP2());
					grid_unpack_with_prp<op,prp_object,device_grid,Memory>::template unpacking<decltype(sub2),prp...>(prRecv_prp,sub2,loc_grid.get(sub_id),ps);
				}
			}
		}
		else
		{
			// Unpack the object
			for ( size_t i = 0 ; i < recv_buffers.size() ; i++ )
			{
				Unpack_stat ps;
				size_t mark_here = ps.getOffset();

				ExtPreAlloc<BMemory<HeapMemory>> mem(recv_buffers.get(i).size(),recv_buffers.get(i));

				// for each external ghost box
				while (ps.getOffset() - mark_here < recv_buffers.get(i).size())
				{
					// Unpack the ghost box global-id

					// Unpack the ghost box global-id

					size_t g_id;
					Unpacker<size_t,BMemory<HeapMemory>>::unpack(mem,g_id,ps);

					size_t pid = dec.ProctoID(recv_proc.get(i).p_id);

					size_t l_id = 0;
					// convert the global id into local id
					auto key = g_id_to_internal_ghost_box.get(pid).find(g_id);
					if (key != g_id_to_internal_ghost_box.get(pid).end()) // FOUND
					{l_id = key->second;}
					else
					{
						// NOT FOUND

						// It must be always found, if not it mean that the processor has no-idea of
						// what is stored and conseguently do not know how to unpack, print a critical error
						// and return

						std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Critical, cannot unpack object, because received data cannot be interpreted\n";

						return;
					}

					// Get the internal ghost box associated with the packed information
					Box<dim,size_t> box = ig_box.get(pid).bid.get(l_id).box;
					size_t sub_id = ig_box.get(pid).bid.get(l_id).sub;
					box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

					// sub-grid where to unpack
					auto sub2 = loc_grid.get(sub_id).getIterator(box.getKP1(),box.getKP2());
					grid_unpack_with_prp<op,prp_object,device_grid,BMemory<HeapMemory>>::template unpacking<decltype(sub2),prp...>(mem,sub2,loc_grid.get(sub_id),ps);
				}
			}
		}
	}

	int find_local_sub(Box<dim, long int> & box_dst, openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext)
	{
		Point<dim,long int> point;
		for (size_t n = 0; n < dim; n++)
		{point.get(n) = (box_dst.getHigh(n) + box_dst.getLow(n))/2;}

		for (size_t j = 0; j < gdb_ext.size(); j++)
		{
			// Local sub-domain
			SpaceBox<dim,long int> sub = gdb_ext.get(j).Dbox;
			sub += gdb_ext.get(j).origin;

			if (sub.isInside(point) == true)
			{
				return j;
			}
		}
		return -1;
	}

public:

	/*! \brief Reconstruct the local grids
	 *
	 * \param m_oGrid_recv Vector of labeled grids to combine into a local grid
	 * \param loc_grid local grids
	 * \param gdb_ext information of the local grids
	 * \param cd_sm Cell-decomposer
	 *
	 */
	inline void grids_reconstruct(openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> & m_oGrid_recv,
			                      openfpm::vector<device_grid> & loc_grid,
								  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
								  CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm)
	{
		// Clear the information of the grid
		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{loc_grid.get(i).clear();}

		for (size_t a = 0; a < m_oGrid_recv.size(); a++)
		{
			for (size_t k = 0; k < m_oGrid_recv.get(a).size(); k++)
			{
				device_grid & g = m_oGrid_recv.get(a).template get<0>(k);

				SpaceBox<dim,long int> b = m_oGrid_recv.get(a).template get<1>(k);

				Point<dim,St> p;
				for (size_t n = 0; n < dim; n++)
				{p.get(n) = g.getGrid().getBox().getHigh(n);}

				Point<dim,St> point;
				for (size_t n = 0; n < dim; n++)
				{point.get(n) = (b.getHigh(n) + b.getLow(n))/2;}

				for (size_t j = 0; j < gdb_ext.size(); j++)
				{
					// Local sub-domain
					SpaceBox<dim,long int> sub = gdb_ext.get(j).Dbox;
					sub += gdb_ext.get(j).origin;

					if (sub.isInside(point) == true)
					{


						grid_key_dx<dim> start = b.getKP1() - grid_key_dx<dim>(gdb_ext.get(j).origin.asArray());
						grid_key_dx<dim> stop = b.getKP2() - grid_key_dx<dim>(gdb_ext.get(j).origin.asArray());

						Box<dim,size_t> box_src;
						Box<dim,size_t> box_dst;

						for(size_t i = 0 ; i < dim ; i++)
						{
							box_dst.setLow(i,start.get(i));
							box_dst.setHigh(i,stop.get(i));
							box_src.setLow(i,0);
							box_src.setHigh(i,stop.get(i)-start.get(i));
						}

						loc_grid.get(j).copy_to(g,box_src,box_dst);
					}
				}
			}
		}

		std::cout << "UNPACKING " << std::endl;

		for (size_t i = 0 ; i < m_oGrid_recv.size() ; i++)
		{
			for (size_t j = 0 ; j < m_oGrid_recv.get(i).size() ; j++)
			{
				m_oGrid_recv.get(i).template get<0>(j).template deviceToHost<0>();
				std::cout << "UNPACKING POINTS: " << m_oGrid_recv.get(i).template get<0>(j).size() << std::endl;
				m_oGrid_recv.get(i).template get<0>(j).template removeCopyToFinalize<0>(v_cl.getgpuContext(), rem_copy_opt::PHASE1);
			}
		}

		for (size_t i = 0 ; i < m_oGrid_recv.size() ; i++)
		{
			for (size_t j = 0 ; j < m_oGrid_recv.get(i).size() ; j++)
			{m_oGrid_recv.get(i).template get<0>(j).template removeCopyToFinalize<0>(v_cl.getgpuContext(), rem_copy_opt::PHASE2);}
		}

		for (size_t i = 0 ; i < m_oGrid_recv.size() ; i++)
		{
			for (size_t j = 0 ; j < m_oGrid_recv.get(i).size() ; j++)
			{m_oGrid_recv.get(i).template get<0>(j).template removeCopyToFinalize<0>(v_cl.getgpuContext(), rem_copy_opt::PHASE3);}
		}
	}

	/*! \brief Label intersection grids for mappings
	 *
	 * \param dec Decomposition
	 * \param loc_grid_old old local grids
	 * \param cd_sm Cell-decomposer
	 * \param gdb_ext information of the local grids
	 * \param gdb_ext_old information of the old local grids
	 * \param gdb_ext_global information of the grids globaly
	 * \param lbl_b label for each grid
	 * \param prc_sz For each processor the number of grids to send to
	 *
	 */
	template<typename lambda_t>
	inline void labelIntersectionGridsProcessor_and_pack(Decomposition & dec,
												CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm,
												openfpm::vector<device_grid> & loc_grid_old,
												openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
												openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old,
												openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global,
												size_t p_id_cur,
												lambda_t f)
	{
		// lbl_bc.clear();
		// lbl_bc.resize(v_cl.getProcessingUnits());

		// for (int i = 0 ; i < lbl_bc.size() ; i++) 
		// {lbl_bc.get(i) = 0;}

		// // count

		// for (size_t i = 0; i < gdb_ext_old.size(); i++)
		// {
		// 	// Local old sub-domain in global coordinates
		// 	SpaceBox<dim,long int> sub_dom = gdb_ext_old.get(i).Dbox;
		// 	sub_dom += gdb_ext_old.get(i).origin;

		// 	for (size_t j = 0; j < gdb_ext_global.size(); j++)
		// 	{
		// 		size_t p_id = 0;

		// 		// Intersection box
		// 		SpaceBox<dim,long int> inte_box;

		// 		// Global new sub-domain in global coordinates
		// 		SpaceBox<dim,long int> sub_dom_new = gdb_ext_global.get(j).Dbox;
		// 		sub_dom_new += gdb_ext_global.get(j).origin;

		// 		bool intersect = false;

		// 		if (sub_dom.isValid() == true && sub_dom_new.isValid() == true)
		// 			intersect = sub_dom.Intersect(sub_dom_new, inte_box);

		// 		if (intersect == true)
		// 		{
		// 			auto inte_box_cont = cd_sm.convertCellUnitsIntoDomainSpace(inte_box);

		// 			// Get processor ID that store intersection box
		// 			Point<dim,St> p;
		// 			for (size_t n = 0; n < dim; n++)
		// 				p.get(n) = (inte_box_cont.getHigh(n) + inte_box_cont.getLow(n))/2;

		// 			p_id = dec.processorID(p);

		// 			lbl_bc.get(p_id) += 1;
		// 		}
		// 	}
		// }

		// // reserve
		// for (int i = 0 ; i < lbl_b.size() ; i++) 
		// {lbl_b.get(i).reserve(lbl_bc.get(i));}


		// Label all the intersection grids with the processor id where they should go

		for (size_t i = 0; i < gdb_ext_old.size(); i++)
		{
			// Local old sub-domain in global coordinates
			SpaceBox<dim,long int> sub_dom = gdb_ext_old.get(i).Dbox;
			sub_dom += gdb_ext_old.get(i).origin;

			for (size_t j = 0; j < gdb_ext_global.size(); j++)
			{
				size_t p_id = 0;

				// Intersection box
				SpaceBox<dim,long int> inte_box;

				// Global new sub-domain in global coordinates
				SpaceBox<dim,long int> sub_dom_new = gdb_ext_global.get(j).Dbox;
				sub_dom_new += gdb_ext_global.get(j).origin;

				bool intersect = false;

				if (sub_dom.isValid() == true && sub_dom_new.isValid() == true)
					intersect = sub_dom.Intersect(sub_dom_new, inte_box);

				if (intersect == true)
				{
					auto inte_box_cont = cd_sm.convertCellUnitsIntoDomainSpace(inte_box);

					// Get processor ID that store intersection box
					Point<dim,St> p;
					for (size_t n = 0; n < dim; n++)
						p.get(n) = (inte_box_cont.getHigh(n) + inte_box_cont.getLow(n))/2;

					p_id = dec.processorID(p);
					if (p_id != p_id_cur)
					{continue;}
//					prc_sz.get(p_id)++;

					// Transform coordinates to local
					auto inte_box_local = inte_box;

					inte_box_local -= gdb_ext_old.get(i).origin;

					// Grid corresponding for gdb_ext_old.get(i) box
					device_grid & gr = loc_grid_old.get(i);

					// Size of the grid to send
					size_t sz[dim];
					for (size_t l = 0; l < dim; l++)
					{
						sz[l] = inte_box_local.getHigh(l) - inte_box_local.getLow(l) + 1;
						//std::cout << "GR_send size on " << l << " dimension: " << sz[l] << std::endl;
					}

					// Grid to send
					//device_grid gr_send(sz);
					//gr_send.setMemory();
					// lbl_b.get(p_id).add();
					// device_grid & gr_send = lbl_b.get(p_id).last().template get<0>();
					// SpaceBox<dim,long int> & box_send = lbl_b.get(p_id).last().template get<1>();
					// gr_send.setMemory();

					// Sub iterator across intersection box inside local grid
					grid_key_dx<dim> start = inte_box_local.getKP1();
					grid_key_dx<dim> stop = inte_box_local.getKP2();

					Box<dim,long int> box_src;
					Box<dim,long int> box_dst;

					for(size_t i = 0 ; i < dim ; i++)
					{
						box_src.setLow(i,start.get(i));
						box_src.setHigh(i,stop.get(i));
						box_dst.setLow(i,inte_box.getLow(i));
						box_dst.setHigh(i,inte_box.getHigh(i));
					}

					f(box_src,box_dst,gr,p_id);
				}
			}
		}

/*		for (size_t i = 0 ; i < loc_grid_old.size() ; i++)
		{
			loc_grid_old.get(i).template removeCopyToFinalize<0>(v_cl.getgpuContext(), rem_copy_opt::PHASE1);
		}

		for (size_t i = 0 ; i < loc_grid_old.size() ; i++)
		{
			loc_grid_old.get(i).template removeCopyToFinalize<0>(v_cl.getgpuContext(), rem_copy_opt::PHASE2);
		}

		for (size_t i = 0 ; i < loc_grid_old.size() ; i++)
		{
			loc_grid_old.get(i).template removeCopyToFinalize<0>(v_cl.getgpuContext(), rem_copy_opt::PHASE3);
		}*/
	}

	/*! \brief Unpack 
	 *
	 *
	 * 
	 */
	template<int ... prp>
	void unpack_buffer_to_local_grid(openfpm::vector<device_grid> & loc_grid,
									 openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
									 ExtPreAlloc<Memory> & send_buffer,
									 size_t sz)
	{
		// unpack local
		Unpack_stat ps;

		while (ps.getOffset() < sz)
		{
			send_buffer.reset();

			Box<dim,long int> box_dst;
			send_buffer.deviceToHost(ps.getOffset(),ps.getOffset()+sizeof(Box<dim,long int>));
			Unpacker<Box<dim,long int>,Memory>::unpack(send_buffer,box_dst,ps);

			int s = find_local_sub(box_dst,gdb_ext);
			if (s == -1)
			{std::cout << __FILE__ << ":" << __LINE__ << " map, error non-local subdomain " << std::endl;}

			// convert box_dst to local
			for (int d = 0 ; d < dim ; d++ )
			{
				box_dst.setLow(d, box_dst.getLow(d) - gdb_ext.get(s).origin.get(d));
				box_dst.setHigh(d, box_dst.getHigh(d) - gdb_ext.get(s).origin.get(d));
			}
			
			loc_grid.get(s).remove(box_dst);
			auto sub2 = loc_grid.get(s).getIterator(box_dst.getKP1(),box_dst.getKP2(),0);
			Unpacker<device_grid,Memory>::template unpack<decltype(sub2),decltype(v_cl.getgpuContext()),prp ...>(send_buffer,sub2,loc_grid.get(s),ps,v_cl.getgpuContext(),NONE_OPT);
		}

		for (int s = 0 ; s < loc_grid.size() ; s++)
		{loc_grid.get(s).template removeAddUnpackFinalize<prp ...>(v_cl.getgpuContext(),0);}
	}

	/*! \brief Moves all the grids that does not belong to the local processor to the respective processor
	 *
	 * This function in general is called if the decomposition change
	 *
	 * \param dec Decomposition
	 * \param cd_sm cell-decomposer
	 * \param loc_grid set of local grids
	 * \param loc_grid_old set of old local grids
	 * \param gdb_ext information of the local grids
	 * \param gdb_ext_old information of the old local grids
	 * \param gdb_ext_global it contain the decomposition at global level
	 *
	 */
	template<int ... prp>
	void map_(Decomposition & dec,
			  CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm,
			  openfpm::vector<device_grid> & loc_grid,
			  openfpm::vector<device_grid> & loc_grid_old,
			  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
			  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old,
			  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global,
			  size_t opt)
	{
		this->opt = opt;

		openfpm::vector<size_t> send_buffer_sizes(v_cl.getProcessingUnits());
		openfpm::vector<Memory> send_buffers_;
		openfpm::vector<ExtPreAlloc<Memory>> send_buffers;
		send_buffers_.resize(v_cl.getProcessingUnits());
		send_buffers.resize(v_cl.getProcessingUnits());

		send_prc_queue.clear();
		send_pointer.clear();
		send_size.clear();

		for (int p_id = 0 ; p_id < v_cl.getProcessingUnits() ; p_id++)
		{
			for (int i = 0 ; i < loc_grid_old.size() ; i++)
			{loc_grid_old.get(i).packReset();}

			auto l = [&](Box<dim,long int> & box_src,
						Box<dim,long int> & box_dst,
						device_grid & gr,
						size_t p_id){
						//gr_send.copy_to(gr,box_src,box_dst);

			
						Packer<SpaceBox<dim,long int>,BMemory<Memory>>::packRequest(box_dst,send_buffer_sizes.get(p_id));

						auto sub_it = gr.getIterator(box_src.getKP1(),box_src.getKP2(),0);
						gr.template packRequest<prp ...>(sub_it,send_buffer_sizes.get(p_id));

						//box_send = inte_box;
			};

			// Contains the processor id of each box (basically where they have to go)
			labelIntersectionGridsProcessor_and_pack(dec,cd_sm,loc_grid_old,gdb_ext,gdb_ext_old,gdb_ext_global,p_id,l);

			for (int i = 0 ; i < loc_grid_old.size(); i++)
			{
				loc_grid_old.get(i).template packCalculate<prp ...>(send_buffer_sizes.get(p_id),v_cl.getgpuContext());
			}

			send_buffers_.get(p_id).resize(send_buffer_sizes.get(p_id));
			send_buffers.get(p_id).setMemory(send_buffer_sizes.get(p_id),send_buffers_.get(p_id));
			send_buffers.get(p_id).incRef();

			// we now pack
			Pack_stat sts;

			auto lp = [&](Box<dim,long int> & box_src,
						Box<dim,long int> & box_dst,
						device_grid & gr,
						size_t p_id){

						size_t offset = send_buffers.get(p_id).getOffsetEnd();
						Packer<Box<dim,long int>,Memory>::pack(send_buffers.get(p_id),box_dst,sts);
						size_t offset2 = send_buffers.get(p_id).getOffsetEnd();

						send_buffers.get(p_id).hostToDevice(offset,offset2);

						auto sub_it = gr.getIterator(box_src.getKP1(),box_src.getKP2(),0);

						Packer<device_grid,Memory>::template pack<decltype(sub_it),prp ...>(send_buffers.get(p_id),gr,sub_it,sts);
			};

			// Contains the processor id of each box (basically where they have to go)
			labelIntersectionGridsProcessor_and_pack(dec,cd_sm,loc_grid_old,gdb_ext,gdb_ext_old,gdb_ext_global,p_id,lp);

			for (int i = 0 ; i < loc_grid_old.size() ; i++)
			{
				loc_grid_old.get(i).template packFinalize<prp ...>(send_buffers.get(p_id),sts,0,false);
			}
		}

		// std::cout << "Local buffer: " << send_buffers.get(v_cl.rank()).size() << std::endl;
		// int sz = send_buffers.get(v_cl.rank()).size();
		//send_buffers.get(v_cl.rank()).reset();

		// // Print all the byte in send_buffers_
		// for (int j = 0 ; j < 16 && j < sz ; j++) {
		// 	std::cout << "Local buffer " << v_cl.rank() << " " << ((long int *)send_buffers.get(v_cl.rank()).getPointer())[j] << " " << &((long int *)send_buffers.get(v_cl.rank()).getPointer())[j] << std::endl;
		// }

		unpack_buffer_to_local_grid<prp ...>(loc_grid,gdb_ext,send_buffers.get(v_cl.rank()),send_buffers.get(v_cl.rank()).size());

		//openfpm::vector<void *> send_pointer;
		//openfpm::vector<int> send_size;
		for (int i = 0 ; i < send_buffers.size() ; i++)
		{
			if (i != v_cl.rank())
			{
				send_pointer.add(send_buffers_.get(i).getDevicePointer());
				send_size.add(send_buffers_.get(i).size());
				send_prc_queue.add(i);
			}
		}

		size_t * send_size_ptr = NULL;
		size_t * send_prc_queue_ptr = NULL;
		void ** send_pointer_ptr = NULL;

		if (send_size.size() != 0)
		{
			send_size_ptr = &send_size.get(0);
			send_pointer_ptr = &send_pointer.get(0);
			send_prc_queue_ptr = &send_prc_queue.get(0);
		}

		recv_buffers.clear();
		recv_proc.clear();

		v_cl.sendrecvMultipleMessagesNBX(send_pointer.size(),send_size_ptr,
											 send_prc_queue_ptr,send_pointer_ptr,
											 receive_dynamic,this);


		for (int i = 0 ; i < recv_buffers.size() ; i++)
		{
			ExtPreAlloc<Memory> prAlloc_;
			prAlloc_.setMemory(recv_buffers.get(i).size(),recv_buffers.get(i));
			unpack_buffer_to_local_grid<prp ...>(loc_grid,gdb_ext,prAlloc_,recv_proc.get(i).size);
		}

		for (int i = 0 ; i < send_buffers.size() ; i++)
		{send_buffers.get(i).decRef();}
	}

	/*! \brief It fill the ghost part of the grids
	 *
	 * \param ig_box internal ghost box
	 * \param eg_box external ghost box
	 * \param loc_ig_box local internal ghost box
	 * \param loc_eg_box local external ghost box
	 * \param gdb_ext local grids information
	 * \param loc_grid set of local grid
	 * \param g_id_to_external_ghost_box index to external ghost box
	 *
	 */
	template<int... prp> void ghost_get_(const openfpm::vector<ip_box_grid<dim>> & ig_box,
									     const openfpm::vector<ep_box_grid<dim>> & eg_box,
										 const openfpm::vector<i_lbox_grid<dim>> & loc_ig_box,
										 const openfpm::vector<e_lbox_grid<dim>> & loc_eg_box,
			                             const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
										 const openfpm::vector<e_box_multi<dim>> & eb_gid_list,
										 bool use_bx_def,
										 openfpm::vector<device_grid> & loc_grid,
										 const grid_sm<dim,void> & ginfo,
										 std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box,
										 size_t opt)
	{
#ifdef PROFILE_SCOREP
		SCOREP_USER_REGION("ghost_get",SCOREP_USER_REGION_TYPE_FUNCTION)
#endif

		// Sending property object
                typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		recv_buffers.clear();
		recv_proc.clear();
		send_prc_queue.clear();
		send_pointer.clear();
		send_size.clear();

		this->opt = opt;

		size_t req = 0;

		// Pack information
		Pack_stat sts;

		// We check if skip labelling is possible in this condition
		for (int i = 0 ; i < loc_grid.size() ; i++)
		{opt &= (loc_grid.get(i).isSkipLabellingPossible())?(int)-1:~SKIP_LABELLING;}

		#ifdef ENABLE_GRID_DIST_ID_PERF_STATS
		timer packing_time;
		packing_time.start();
		#endif

		if (!(opt & SKIP_LABELLING))
		{
			// first we initialize the pack buffer on all internal grids

			for (size_t i = 0 ; i < loc_grid.size() ; i++)
			{loc_grid.get(i).packReset();}

			// Calculating the size to pack all the data to send
			for ( size_t i = 0 ; i < ig_box.size() ; i++ )
			{
				// for each ghost box
				for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
				{
					// And linked sub-domain
					size_t sub_id = ig_box.get(i).bid.get(j).sub;
					// Internal ghost box
					Box<dim,long int> g_ig_box = ig_box.get(i).bid.get(j).box;

					if (g_ig_box.isValid() == false)
					{continue;}

					g_ig_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

					// Pack a size_t for the internal ghost id
					Packer<size_t,Memory>::packRequest(req);
					// Create a sub grid iterator spanning the internal ghost layer
					auto sub_it = loc_grid.get(sub_id).getIterator(g_ig_box.getKP1(),g_ig_box.getKP2(),false);

					// get the size to pack
					Packer<device_grid,Memory>::template packRequest<decltype(sub_it),prp...>(loc_grid.get(sub_id),sub_it,req);
				}
			}

			// Finalize calculation
			for (size_t i = 0 ; i < loc_grid.size() ; i++)
			{loc_grid.get(i).template packCalculate<prp ...>(req,v_cl.getgpuContext());}

			// resize the property buffer memory
			g_send_prp_mem.resize(req);

			// Create an object of preallocated memory for properties
			ExtPreAlloc<Memory> & prAlloc_prp = *(new ExtPreAlloc<Memory>(req,g_send_prp_mem));
			// Necessary. We do not want this memory to be destroyed untill is going out of scope
			// P.S. Packer shaoe this memory with data-structures and data structures if they see the
			// reference counter to zero they destriy this memory
			prAlloc_prp.incRef();

			pointers.clear();
			pointers2.clear();

			// Pack the information for each processor and send it
			for ( size_t i = 0 ; i < ig_box.size() ; i++ )
			{

				sts.mark();

				void * pointer;

				if (opt & RUN_ON_DEVICE)
				{pointer = prAlloc_prp.getDevicePointerEnd();}
				else
				{pointer = prAlloc_prp.getPointerEnd();}

				// for each ghost box
				for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
				{
					// we pack only if it is valid
					if (ig_box.get(i).bid.get(j).box.isValid() == false)
						continue;

					// And linked sub-domain
					size_t sub_id = ig_box.get(i).bid.get(j).sub;
					// Internal ghost box
					Box<dim,size_t> g_ig_box = ig_box.get(i).bid.get(j).box;
					g_ig_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();
					// Ghost box global id
					size_t g_id = ig_box.get(i).bid.get(j).g_id;

					// Pack a size_t for the internal ghost id
					Packer<size_t,Memory>::pack(prAlloc_prp,g_id,sts);
					prAlloc_prp.hostToDevice(prAlloc_prp.getOffset(),prAlloc_prp.getOffsetEnd());
					// Create a sub grid iterator spanning the internal ghost layer
					auto sub_it = loc_grid.get(sub_id).getIterator(g_ig_box.getKP1(),g_ig_box.getKP2(),false);
					// and pack the internal ghost grid
					Packer<device_grid,Memory>::template pack<decltype(sub_it),prp...>(prAlloc_prp,loc_grid.get(sub_id),sub_it,sts);
				}
				// send the request

				void * pointer2;

				if (opt & RUN_ON_DEVICE)
				{pointer2 = prAlloc_prp.getDevicePointerEnd();}
				else
				{pointer2 = prAlloc_prp.getPointerEnd();}

				pointers.add(pointer);
				pointers2.add(pointer2);
			}

			for (size_t i = 0 ; i < loc_grid.size() ; i++)
			{
				rem_copy_opt opt_ = rem_copy_opt::NONE_OPT;
				if (opt & SKIP_LABELLING)
				{opt_ = rem_copy_opt::KEEP_GEOMETRY;}

				loc_grid.get(i).template packFinalize<prp ...>(prAlloc_prp,sts,opt_,true);
			}

			prAlloc_prp.decRef();
			delete &prAlloc_prp;
		}
		else
		{
			req = g_send_prp_mem.size();

			// Create an object of preallocated memory for properties
			ExtPreAlloc<Memory> & prAlloc_prp = *(new ExtPreAlloc<Memory>(req,g_send_prp_mem));
			prAlloc_prp.incRef();

			for (size_t i = 0 ; i < loc_grid.size() ; i++)
			{
				rem_copy_opt opt_ = rem_copy_opt::NONE_OPT;
				if (opt & SKIP_LABELLING)
				{opt_ = rem_copy_opt::KEEP_GEOMETRY;}

				loc_grid.get(i).template packFinalize<prp ...>(prAlloc_prp,sts,opt_,true);
			}

			prAlloc_prp.decRef();
			delete &prAlloc_prp;
		}

		#ifdef ENABLE_GRID_DIST_ID_PERF_STATS
		packing_time.stop();
		tot_pack += packing_time.getwct();
		timer sendrecv_time;
		sendrecv_time.start();
		#endif

		for ( size_t i = 0 ; i < ig_box.size() ; i++ )
		{
			// This function send (or queue for sending) the information
			send_or_queue(ig_box.get(i).prc,(char *)pointers.get(i),(char *)pointers2.get(i));
		}

		// Calculate the total information to receive from each processors
		std::vector<size_t> prp_recv;

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prRecv_prp = *(new ExtPreAlloc<Memory>(g_recv_prp_mem.size(),g_recv_prp_mem));
		prRecv_prp.incRef();

		// Before wait for the communication to complete we sync the local ghost
		// in order to overlap with communication

		queue_recv_data_get<prp_object>(eg_box,prp_recv,prRecv_prp);

		#ifdef ENABLE_GRID_DIST_ID_PERF_STATS
		sendrecv_time.stop();
		tot_sendrecv += sendrecv_time.getwct();
		timer merge_loc_time;
		merge_loc_time.start();
		#endif

		ghost_get_local<prp...>(loc_ig_box,loc_eg_box,gdb_ext,loc_grid,g_id_to_external_ghost_box,ginfo,use_bx_def,opt);

		#ifdef ENABLE_GRID_DIST_ID_PERF_STATS
		merge_loc_time.stop();
		tot_loc_merge += merge_loc_time.getwct();
		timer merge_time;
		merge_time.start();
		#endif

		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{loc_grid.get(i).removeAddUnpackReset();}

		merge_received_data_get<prp ...>(loc_grid,eg_box,prp_recv,prRecv_prp,g_id_to_external_ghost_box,eb_gid_list,opt);

		rem_copy_opt opt_ = rem_copy_opt::NONE_OPT;
		if (opt & SKIP_LABELLING)
		{opt_ = rem_copy_opt::KEEP_GEOMETRY;}

		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{loc_grid.get(i).template removeAddUnpackFinalize<prp ...>(v_cl.getgpuContext(),opt_);}

		#ifdef ENABLE_GRID_DIST_ID_PERF_STATS
		merge_time.stop();
		tot_merge += merge_time.getwct();
		#endif

		prRecv_prp.decRef();
		delete &prRecv_prp;
	}

	/*! \brief It merge the information in the ghost with the
	 *         real information
	 *
	 * \tparam op merge operation
	 *
	 * \param ig_box internal ghost box
	 * \param eg_box external ghost box
	 * \param loc_ig_box local internal ghost box
	 * \param loc_eg_box local external ghost box
	 * \param gdb_ext local grids information
	 * \param loc_grid set of local grid
	 * \param g_id_to_internal_ghost_box index to internal ghost box
	 *
	 */
	template<template<typename,typename> class op,int... prp>
	void ghost_put_(Decomposition & dec,
			        const openfpm::vector<ip_box_grid<dim>> & ig_box,
					const openfpm::vector<ep_box_grid<dim>> & eg_box,
					const openfpm::vector<i_lbox_grid<dim>> & loc_ig_box,
					const openfpm::vector<e_lbox_grid<dim>> & loc_eg_box,
			        const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
					openfpm::vector<device_grid> & loc_grid,
					openfpm::vector<std::unordered_map<size_t,size_t>> & g_id_to_internal_ghost_box)
	{
		// Sending property object
		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		recv_buffers.clear();
		recv_proc.clear();
		send_prc_queue.clear();
		send_pointer.clear();
		send_size.clear();

		size_t req = 0;

		// Create a packing request vector
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			// for each ghost box
			for (size_t j = 0 ; j < eg_box.get(i).bid.size() ; j++)
			{
				// And linked sub-domain
				size_t sub_id = eg_box.get(i).bid.get(j).sub;
				// Internal ghost box
				Box<dim,long int> g_eg_box = eg_box.get(i).bid.get(j).g_e_box;

				if (g_eg_box.isValid() == false)
					continue;

				g_eg_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

				// Pack a size_t for the internal ghost id
				Packer<size_t,HeapMemory>::packRequest(req);

				// Create a sub grid iterator spanning the external ghost layer
				auto sub_it = loc_grid.get(sub_id).getIterator(g_eg_box.getKP1(),g_eg_box.getKP2());

				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template packRequest<decltype(sub_it),prp...>(loc_grid.get(sub_id),sub_it,req);
			}
		}

		// resize the property buffer memory
		g_send_prp_mem.resize(req);

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prAlloc_prp = *(new ExtPreAlloc<Memory>(req,g_send_prp_mem));

		prAlloc_prp.incRef();

		// Pack information
		Pack_stat sts;

		// Pack the information for each processor and send it
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{

			sts.mark();
			void * pointer = prAlloc_prp.getPointerEnd();

			// for each ghost box
			for (size_t j = 0 ; j < eg_box.get(i).bid.size() ; j++)
			{
				// we pack only if it is valid
				if (eg_box.get(i).bid.get(j).g_e_box.isValid() == false)
					continue;

				// And linked sub-domain
				size_t sub_id = eg_box.get(i).bid.get(j).sub;
				// Internal ghost box
				Box<dim,size_t> g_eg_box = eg_box.get(i).bid.get(j).g_e_box;
				g_eg_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();
				// Ghost box global id
				size_t g_id = eg_box.get(i).bid.get(j).g_id;

				// Pack a size_t for the internal ghost id
				Packer<size_t,HeapMemory>::pack(prAlloc_prp,g_id,sts);
				// Create a sub grid iterator spanning the external ghost layer
				auto sub_it = loc_grid.get(sub_id).getIterator(g_eg_box.getKP1(),g_eg_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template pack<decltype(sub_it),prp...>(prAlloc_prp,loc_grid.get(sub_id),sub_it,sts);
			}
			// send the request

			void * pointer2 = prAlloc_prp.getPointerEnd();

			// This function send (or queue for sending) the information
			send_or_queue(ig_box.get(i).prc,(char *)pointer,(char *)pointer2);
		}

		// Calculate the total information to receive from each processors
		std::vector<size_t> prp_recv;

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prRecv_prp = *(new ExtPreAlloc<Memory>(0,g_recv_prp_mem));
		prRecv_prp.incRef();

		queue_recv_data_put<prp_object>(ig_box,prp_recv,prRecv_prp);

		// Before wait for the communication to complete we sync the local ghost
		// in order to overlap with communication

		ghost_put_local<op,prp...>(loc_ig_box,loc_eg_box,gdb_ext,loc_grid,g_id_to_internal_ghost_box);

		merge_received_data_put<op,prp ...>(dec,loc_grid,ig_box,prp_recv,prRecv_prp,gdb_ext,g_id_to_internal_ghost_box);

		prRecv_prp.decRef();
		prAlloc_prp.decRef();
		delete &prAlloc_prp;
		delete &prRecv_prp;
	}

	/*! \brief Constructor
	 *
	 *
	 */
	grid_dist_id_comm()
	:v_cl(create_vcluster<Memory>())
	{

	}

	/*! \brief Copy constructor
	 *
	 * It does not really copy. This structure it suppose to store only
	 * temporal data
	 *
	 */
	grid_dist_id_comm(const grid_dist_id_comm<dim,St,T,Decomposition,Memory,device_grid> & gc)
	:v_cl(gc.v_cl)
	{

	}
};


#endif /* SRC_GRID_GRID_DIST_ID_COMM_HPP_ */
