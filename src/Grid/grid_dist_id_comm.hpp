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
	template<template<typename,typename> class op, int ... prp> static void call_unpack(ExtPreAlloc<Memory> & recv_buf, grid_key_dx_iterator_sub<device_grid::dims> & sub2, device_grid & gd, Unpack_stat & ps)
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
	template<template<typename,typename> class op, unsigned int ... prp> static void call_unpack(ExtPreAlloc<Memory> & recv_buf, grid_key_dx_iterator_sub<device_grid::dims> & sub2, device_grid & gd, Unpack_stat & ps)
	{
		PtrMemory * ptr1;

		size_t sz[device_grid::dims];

		for (size_t i = 0 ; i < device_grid::dims ; i++)
			sz[i] = sub2.getStop().get(i) - sub2.getStart().get(i) + 1;

		size_t tot = 1;

		for (size_t i = 0 ; i < device_grid::dims ; i++)
			tot *= sz[i];

		tot *= sizeof(T);

#ifdef SE_CLASS1

		if (ps.getOffset() + tot > recv_buf.size())
			std::cerr << __FILE__ << ":" << __LINE__ << " Error: overflow in the receiving buffer for ghost_put" << std::endl;

#endif

		// add the received particles to the vector
		ptr1 = new PtrMemory(((char *)recv_buf.getPointerBase()+ps.getOffset()),tot);

		// create vector representation to a piece of memory already allocated
		grid_cpu<device_grid::dims,T,PtrMemory,typename memory_traits_lin<T>::type> gs;

		gs.setMemory(*ptr1);

		// resize with the number of elements
		gs.resize(sz);

		// Merge the information

		auto it_src = gs.getIterator();

		while (sub2.isNext())
		{
			object_s_di_op<op,
			            decltype(gs.get_o(it_src.get())),
						decltype(gd.get_o(sub2.get())),
						OBJ_ENCAP,prp...>
			(gs.get_o(it_src.get()),
			 gd.get_o(sub2.get()));

			++sub2;
			++it_src;
		}

		ps.addOffset(tot);
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
	template<template<typename,typename> class op, typename T> inline static void call_unpack(ExtPreAlloc<Memory> & recv_buf, grid_key_dx_iterator_sub<device_grid::dims> & sub2, device_grid & dg, Unpack_stat & ps)
	{
		const bool result = has_pack_gen<typename T::type>::value == false;

		grid_unpack_selector_with_prp<result,T,device_grid,Memory>::template call_unpack<op,prp...>(recv_buf,sub2,dg,ps);
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
	template<unsigned int ... prp> static void unpacking(ExtPreAlloc<Memory> & recv_buf, grid_key_dx_iterator_sub<device_grid::dims> & sub2, device_grid & dg, Unpack_stat & ps)
	{
		typedef index_tuple<prp...> ind_prop_to_pack;
		grid_call_serialize_variadic<device_grid,Memory,ind_prop_to_pack>::template call_unpack<op,T>(recv_buf, sub2, dg, ps);
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
	Vcluster & v_cl;

	//! Maps the processor id with the communication request into map procedure
	openfpm::vector<size_t> p_map_req;

	//! Stores the list of processors that communicate with us (local processor)
	openfpm::vector<size_t> prc_recv_map;

	//! Stores the size of the elements added for each processor that communicate with us (local processor)
	openfpm::vector<size_t> recv_sz_map;

	//! For each near processor, outgoing intersection grid
	//! \warning m_oGrid is assumed to be an ordered list
	//! first id is grid
	//! second id is the processor id
	openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> m_oGrid;

	//! Memory for the ghost sending buffer
	Memory g_send_prp_mem;

	//! Memory for the ghost sending buffer
	Memory g_recv_prp_mem;


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
											  std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box)
	{
		grid_key_dx<dim> cnt[1];
		cnt[0].zero();

		//! For all the sub-domains
		for (size_t i = 0 ; i < loc_ig_box.size() ; i++)
		{
			//! For all the internal ghost boxes of each sub-domain
			for (size_t j = 0 ; j < loc_ig_box.get(i).bid.size() ; j++)
			{
				Box<dim,size_t> bx_src = loc_ig_box.get(i).bid.get(j).box;
				// convert into local
				bx_src -= gdb_ext.get(i).origin;

				// sub domain connected with external box
				size_t sub_id_dst = loc_ig_box.get(i).bid.get(j).sub;

				// local internal ghost box connected
				size_t k = loc_ig_box.get(i).bid.get(j).k;

				Box<dim,size_t> bx_dst = loc_eg_box.get(sub_id_dst).bid.get(k).box;

				// convert into local
				bx_dst -= gdb_ext.get(sub_id_dst).origin;

				// create 2 sub grid iterator

				if (bx_dst.isValid() == false)
					continue;

				const auto & gs = loc_grid.get(i);
				auto & gd = loc_grid.get(sub_id_dst);

#ifdef SE_CLASS1

				if (loc_eg_box.get(sub_id_dst).bid.get(k).sub != i)
				{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination are not correctly linked" << "\n";}

				if (bx_src.getVolumeKey() != bx_dst.getVolumeKey())
				{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination does not match in size" << "\n";}

				auto bxs = gs.getGrid().getBoxKey();
				auto bxd = gd.getGrid().getBoxKey();

				if (bxs.isContained(bx_src) == false)
				{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the source box is out of bound of the local grid" << "\n";}

				if (bxd.isContained(bx_dst) == false)
				{std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " the destination box is out of bound of the local grid" << "\n";}

#endif

				typedef typename std::remove_reference<decltype(gd)>::type grid_cp;
				typedef typename std::remove_reference<decltype(loc_grid.get(i).getGrid())>::type grid_info_cp;

				copy_grid_fast<!is_contiguos<prp...>::type::value || has_pack_gen<typename device_grid::value_type>::value,
							   dim,
							   grid_cp,
							   grid_info_cp>::copy(loc_grid.get(i).getGrid(),
						       loc_grid.get(sub_id_dst).getGrid(),
							   bx_src,
							   bx_dst,
							   gs,gd,cnt);
			}
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

				Box<dim,size_t> bx_src = loc_eg_box.get(i).bid.get(j).box;
				// convert into local
				bx_src -= gdb_ext.get(i).origin;

				// sub domain connected with external box
				size_t sub_id_dst = loc_eg_box.get(i).bid.get(j).sub;

				// local external ghost box connected
				size_t k = loc_eg_box.get(i).bid.get(j).k;

				Box<dim,size_t> bx_dst = loc_ig_box.get(sub_id_dst).bid.get(k).box;

				// convert into local
				bx_dst -= gdb_ext.get(sub_id_dst).origin;

				// create 2 sub grid iterator

				if (bx_dst.isValid() == false)
					continue;

				grid_key_dx_iterator_sub<dim> sub_src(loc_grid.get(i).getGrid(),bx_src.getKP1(),bx_src.getKP2());
				grid_key_dx_iterator_sub<dim> sub_dst(loc_grid.get(sub_id_dst).getGrid(),bx_dst.getKP1(),bx_dst.getKP2());

#ifdef SE_CLASS1

				if (loc_ig_box.get(sub_id_dst).bid.get(k).sub != i)
					std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination are not correctly linked" << "\n";

				if (sub_src.getVolume() != sub_dst.getVolume())
					std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination does not match in size" << "\n";

#endif

				const auto & gs = loc_grid.get(i);
				auto & gd = loc_grid.get(sub_id_dst);

				while (sub_src.isNext())
				{
					// write the object in the last element
					object_s_di_op<op,decltype(gs.get_o(sub_src.get())),decltype(gd.get_o(sub_dst.get())),OBJ_ENCAP,prp...>(gs.get_o(sub_src.get()),gd.get_o(sub_dst.get()));

					++sub_src;
					++sub_dst;
				}
			}
		}
	}

	/*! \brief this function create send and receive asynchronously to receive ghosts part
	 *
	 * \param ig_box internal ghost box
	 * \param eg_box external ghost box
	 *
	 */
	template<int... prp>
	void send_and_receive_ghost(ExtPreAlloc<Memory> ** prAlloc_prp,
								ExtPreAlloc<Memory> ** prRecv_prp,
								const openfpm::vector<ip_box_grid<dim>> & ig_box,
								const openfpm::vector<ep_box_grid<dim>> & eg_box,
								const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
								openfpm::vector<device_grid> & loc_grid,
								size_t & req)
	{
		// Sending property object
		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		// Create a packing request vector
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
					continue;

				g_ig_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

				// Pack a size_t for the internal ghost id
				Packer<size_t,HeapMemory>::packRequest(req);

				// Create a sub grid iterator spanning the internal ghost layer

				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_ig_box.getKP1(),g_ig_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template packRequest<prp...>(loc_grid.get(sub_id),sub_it,req);
			}
		}

		// resize the property buffer memory
		g_send_prp_mem.resize(req);

		// Create an object of preallocated memory for properties
		(*prAlloc_prp) = new ExtPreAlloc<Memory>(req,g_send_prp_mem);
		(*prAlloc_prp)->incRef();

		// Pack information
		Pack_stat sts;

		// Pack the information for each processor and send it
		for ( size_t i = 0 ; i < ig_box.size() ; i++ )
		{

			sts.mark();
			void * pointer = (*prAlloc_prp)->getPointerEnd();

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
				Packer<size_t,HeapMemory>::pack(**prAlloc_prp,g_id,sts);

				// Create a sub grid iterator spanning the internal ghost layer
				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_ig_box.getKP1(),g_ig_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template pack<prp...>(**prAlloc_prp,loc_grid.get(sub_id),sub_it,sts);
			}
			// send the request

			void * pointer2 = (*prAlloc_prp)->getPointerEnd();

			v_cl.send(ig_box.get(i).prc,0,pointer,(char *)pointer2 - (char *)pointer);
		}

		// Calculate the total information to receive from each processors
		std::vector<size_t> prp_recv;

		//! Receive the information from each processors
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			prp_recv.push_back(0);

			// for each external ghost box
			for (size_t j = 0 ; j < eg_box.get(i).bid.size() ; j++)
			{
				// External ghost box
				Box<dim,size_t> g_eg_box = eg_box.get(i).bid.get(j).g_e_box;
				prp_recv[prp_recv.size()-1] += g_eg_box.getVolumeKey() * sizeof(prp_object) + sizeof(size_t);
			}
		}

		size_t tot_recv = ExtPreAlloc<Memory>::calculateMem(prp_recv);

		//! Resize the receiving buffer
		g_recv_prp_mem.resize(tot_recv);

		// Create an object of preallocated memory for properties
		(*prRecv_prp) = new ExtPreAlloc<Memory>(tot_recv,g_recv_prp_mem);
		(*prRecv_prp)->incRef();

		// queue the receives
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			(*prRecv_prp)->allocate(prp_recv[i]);
			v_cl.recv(eg_box.get(i).prc,0,(*prRecv_prp)->getPointer(),prp_recv[i]);
		}
	}

	/*! \brief Process the received data
	 *
	 * \param eg_box external ghost box
	 *
	 */
	template<int... prp>
	void process_received(ExtPreAlloc<Memory> * prRecv_prp,
			   	   	   	  const openfpm::vector<ep_box_grid<dim>> & eg_box,
						  openfpm::vector<device_grid> & loc_grid,
						  std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box)
	{
		Unpack_stat ps;

		// Unpack the object
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			// for each external ghost box
			for (size_t j = 0 ; j < eg_box.get(i).bid.size() ; j++)
			{
				// Unpack the ghost box global-id

				size_t g_id;
				Unpacker<size_t,HeapMemory>::unpack(*prRecv_prp,g_id,ps);

				size_t l_id = 0;
				// convert the global id into local id
				auto key = g_id_to_external_ghost_box.find(g_id);
				if (key != g_id_to_external_ghost_box.end()) // FOUND
					l_id = key->second;
				else
				{
					// NOT FOUND

					// It must be always found, if not it mean that the processor has no-idea of
					// what is stored and conseguently do not know how to unpack, print a critical error
					// and return

					std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " Critical, cannot unpack object, because received data cannot be interpreted\n";

					return;
				}

				// Get the external ghost box associated with the packed information
				Box<dim,size_t> box = eg_box.get(i).bid.get(l_id).l_e_box;
				size_t sub_id = eg_box.get(i).bid.get(l_id).sub;

				// sub-grid where to unpack
				grid_key_dx_iterator_sub<dim> sub2(loc_grid.get(sub_id).getGrid(),box.getKP1(),box.getKP2());

				// Unpack
				Unpacker<device_grid,HeapMemory>::template unpack<prp...>(*prRecv_prp,sub2,loc_grid.get(sub_id),ps);
			}
		}
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
		size_t count2 = 0;
		for (size_t a = 0; a < m_oGrid_recv.size(); a++)
		{
			for (size_t k = 0; k < m_oGrid_recv.get(a).size(); k++)
			{
				device_grid g = m_oGrid_recv.get(a).template get<0>(k);

				size_t count = 0;


				auto it = g.getIterator();

				while (it.isNext())
				{
					//auto key = it.get();

					//if (g.template get<0>(key) != 1)
						//std::cout << "WRONG???????" << std::endl;

					++it;
					count++;
				}

				SpaceBox<dim,long int> b = m_oGrid_recv.get(a).template get<1>(k);

				Point<dim,St> p;
				for (size_t n = 0; n < dim; n++)
					p.get(n) = g.getGrid().getBox().getHigh(n);

				//std::cout << "G after resize: (" << g.getGrid().getBox().getLow(0) << "; " << g.getGrid().getBox().getLow(1) << "); (" << g.getGrid().getBox().getHigh(0) << "; " << g.getGrid().getBox().getHigh(1) << ")" << std::endl;

				Point<dim,St> point;
				for (size_t n = 0; n < dim; n++)
					point.get(n) = (b.getHigh(n) + b.getLow(n))/2;

				for (size_t j = 0; j < gdb_ext.size(); j++)
				{
					// Local sub-domain
					SpaceBox<dim,long int> sub = gdb_ext.get(j).Dbox;
					sub += gdb_ext.get(j).origin;

					if (sub.isInside(point) == true)
					{
						grid_key_dx<dim> start = b.getKP1() - grid_key_dx<dim>(gdb_ext.get(j).origin.asArray());
						grid_key_dx<dim> stop = b.getKP2() - grid_key_dx<dim>(gdb_ext.get(j).origin.asArray());

						std::string start2 = start.to_string();
						std::string stop2 = stop.to_string();

						auto it = loc_grid.get(j).getSubIterator(start,stop);

						// Copy selected elements into a local grid
						while (it.isNext())
						{
							auto key = it.get();
							std::string str = key.to_string();
							grid_key_dx<dim> key2 = key - start;

							//std::cout << "Key: " << str << std::endl;
							loc_grid.get(j).get_o(key) = g.get_o(key2);
							count2++;

							++it;
						}
					}
				}
			}
		}
		//std::cout << "Count after: " << count2 << std::endl;
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
	inline void labelIntersectionGridsProcessor(Decomposition & dec,
												CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm,
												openfpm::vector<device_grid> & loc_grid_old,
												openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
												openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old,
												openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global,
												openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> & lbl_b,
												openfpm::vector<size_t> & prc_sz)
	{
		// resize the label buffer
		lbl_b.resize(v_cl.getProcessingUnits());

		size_t count2 = 0;

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
					count2++;
					auto inte_box_cont = cd_sm.convertCellUnitsIntoDomainSpace(inte_box);

					// Get processor ID that store intersection box
					Point<dim,St> p;
					for (size_t n = 0; n < dim; n++)
						p.get(n) = (inte_box_cont.getHigh(n) + inte_box_cont.getLow(n))/2;

					p_id = dec.processorID(p);
					prc_sz.get(p_id)++;

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
					device_grid gr_send(sz);
					gr_send.setMemory();

					// Sub iterator across intersection box inside local grid
					grid_key_dx<dim> start = inte_box_local.getKP1();
					grid_key_dx<dim> stop = inte_box_local.getKP2();

					Point<dim,St> p1;
					for (size_t n = 0; n < dim; n++)
						p1.get(n) = gr_send.getGrid().getBox().getLow(n);

					Point<dim,St> p2;
					for (size_t n = 0; n < dim; n++)
						p2.get(n) = gr_send.getGrid().getBox().getHigh(n);

					std::string start2 = start.to_string();
					std::string stop2 = stop.to_string();

					auto it = gr.getSubIterator(start,stop);

					// Copy selected elements into a new sub-grid
					while (it.isNext())
					{
						auto key = it.get();
						grid_key_dx<dim> key2 = key - start;
						std::string str = key.to_string();

						gr_send.get_o(key2) = gr.get_o(key);

						++it;
					}

					aggregate<device_grid,SpaceBox<dim,long int>> aggr;

					aggr.template get<0>() = gr_send;
					aggr.template get<1>() = inte_box;

					// Add to the labeling vector
					lbl_b.get(p_id).add(aggr);
				}
			}
		}
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
	void map_(Decomposition & dec,
			  CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm,
			  openfpm::vector<device_grid> & loc_grid,
			  openfpm::vector<device_grid> & loc_grid_old,
			  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
			  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_old,
			  openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext_global)
	{
		// Processor communication size
		openfpm::vector<size_t> prc_sz(v_cl.getProcessingUnits());

		// Contains the processor id of each box (basically where they have to go)
		labelIntersectionGridsProcessor(dec,cd_sm,loc_grid_old,gdb_ext,gdb_ext_old,gdb_ext_global,m_oGrid,prc_sz);

		// Calculate the sending buffer size for each processor, put this information in
		// a contiguous buffer
		p_map_req.resize(v_cl.getProcessingUnits());

		// Vector of number of sending grids for each involved processor
		openfpm::vector<size_t> prc_sz_r;
		// Vector of ranks of involved processors
		openfpm::vector<size_t> prc_r;

		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (m_oGrid.get(i).size() != 0)
			{
				p_map_req.get(i) = prc_r.size();
				prc_r.add(i);
				prc_sz_r.add(m_oGrid.get(i).size());
			}
		}

		decltype(m_oGrid) m_oGrid_new;
		for (size_t i = 0; i < v_cl.getProcessingUnits(); i++)
		{
			if (m_oGrid.get(i).size() != 0)
				m_oGrid_new.add(m_oGrid.get(i));
		}

		// Vector for receiving of intersection grids
		openfpm::vector<openfpm::vector<aggregate<device_grid,SpaceBox<dim,long int>>>> m_oGrid_recv;

		// Send and recieve intersection grids
		v_cl.SSendRecv(m_oGrid_new,m_oGrid_recv,prc_r,prc_recv_map,recv_sz_map);

		// Reconstruct the new local grids
		grids_reconstruct(m_oGrid_recv,loc_grid,gdb_ext,cd_sm);
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
										 openfpm::vector<device_grid> & loc_grid,
										 std::unordered_map<size_t,size_t> & g_id_to_external_ghost_box)
	{
		size_t req = 0;

		ExtPreAlloc<Memory> * prRecv_prp = NULL;
		ExtPreAlloc<Memory> * prAlloc_prp = NULL;

		if (v_cl.getProcessingUnits() != 1)
		{send_and_receive_ghost<prp...>(&prAlloc_prp,&prRecv_prp, ig_box,eg_box,gdb_ext,loc_grid,req);}

		// Before wait for the communication to complete we sync the local ghost
		// in order to overlap with communication

		ghost_get_local<prp...>(loc_ig_box,loc_eg_box,gdb_ext,loc_grid,g_id_to_external_ghost_box);

		// wait to receive communication
		v_cl.execute();

		if (v_cl.getProcessingUnits() != 1)
		{process_received<prp...>(prRecv_prp,eg_box,loc_grid,g_id_to_external_ghost_box);}
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
	void ghost_put_(const openfpm::vector<ip_box_grid<dim>> & ig_box,
					const openfpm::vector<ep_box_grid<dim>> & eg_box,
					const openfpm::vector<i_lbox_grid<dim>> & loc_ig_box,
					const openfpm::vector<e_lbox_grid<dim>> & loc_eg_box,
			        const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
					openfpm::vector<device_grid> & loc_grid,
					openfpm::vector<std::unordered_map<size_t,size_t>> & g_id_to_internal_ghost_box)
	{
		// Sending property object
		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

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

				// Create a sub grid iterator spanning the internal ghost layer
				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_eg_box.getKP1(),g_eg_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template packRequest<prp...>(loc_grid.get(sub_id),sub_it,req);
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

				// Create a sub grid iterator spanning the internal ghost layer
				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_eg_box.getKP1(),g_eg_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template pack<prp...>(prAlloc_prp,loc_grid.get(sub_id),sub_it,sts);
			}
			// send the request

			void * pointer2 = prAlloc_prp.getPointerEnd();

			v_cl.send(ig_box.get(i).prc,0,pointer,(char *)pointer2 - (char *)pointer);
		}

		// Calculate the total information to receive from each processors
		std::vector<size_t> prp_recv;

		//! Receive the information from each processors
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

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prRecv_prp = *(new ExtPreAlloc<Memory>(tot_recv,g_recv_prp_mem));
		prRecv_prp.incRef();

		// queue the receives
		for ( size_t i = 0 ; i < ig_box.size() ; i++ )
		{
			prRecv_prp.allocate(prp_recv[i]);
			v_cl.recv(ig_box.get(i).prc,0,prRecv_prp.getPointer(),prp_recv[i]);
		}

		// Before wait for the communication to complete we sync the local ghost
		// in order to overlap with communication

		ghost_put_local<op,prp...>(loc_ig_box,loc_eg_box,gdb_ext,loc_grid,g_id_to_internal_ghost_box);

		// wait to receive communication
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
					l_id = key->second;
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
				grid_key_dx_iterator_sub<dim> sub2(loc_grid.get(sub_id).getGrid(),box.getKP1(),box.getKP2());

				grid_unpack_with_prp<op,prp_object,device_grid,Memory>::template unpacking<prp...>(prRecv_prp,sub2,loc_grid.get(sub_id),ps);
			}
		}
	}

	/*! \brief Constructor
	 *
	 *
	 */
	grid_dist_id_comm()
	:v_cl(create_vcluster())
	{

	}
};


#endif /* SRC_GRID_GRID_DIST_ID_COMM_HPP_ */
