#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include <unordered_map>
#include "Grid/map_grid.hpp"
#include "VCluster.hpp"
#include "Space/SpaceBox.hpp"
#include "util/mathutil.hpp"
#include "grid_dist_id_iterator.hpp"
#include "grid_dist_key.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "VTKWriter.hpp"
#include "Packer.hpp"
#include "Unpacker.hpp"

#define SUB_UNIT_FACTOR 64


/*! \brief This is a distributed grid
 *
 * Implementation of a distributed grid with decomposition on the ids.
 * A distributed grid is a grid distributed across processors.
 * The decomposition is performed on the ids of the elements
 *
 *
 * \param dim Dimensionality of the grid
 * \param St Type of space where the grid is living
 * \param T object the grid is storing
 * \param Decomposition Class that decompose the grid for example CartDecomposition
 * \param Mem Is the allocator
 * \param device type of base structure is going to store the data
 *
 */

template<unsigned int dim, typename St, typename T, typename Decomposition,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_id
{
	// Domain
	Box<dim,St> domain;

	// Ghost expansion
	Ghost<dim,St> ghost;

	//! Local grids
	Vcluster_object_array<device_grid> loc_grid;

	//! Space Decomposition
	Decomposition dec;

	//! Extension of each grid: Domain and ghost + domain
	openfpm::vector<GBoxes<device_grid::dims>> gdb_ext;

	//! Size of the grid on each dimension
	size_t g_sz[dim];

	//! Structure that divide the space into cells
	CellDecomposer_sm<dim,St> cd_sm;

	//! Communicator class
	Vcluster & v_cl;

	//! It map a global ghost id (g_id) to the external ghost box information
	std::unordered_map<size_t,size_t> g_id_to_external_ghost_box;

	/*! \brief Get the grid size
	 *
	 * Given a domain, the resolution of the grid on it and another spaceBox contained in the domain
	 * it give the size on all directions of the local grid
	 *
	 * \param sp SpaceBox enclosing the local grid
	 * \param domain Space box enclosing the physical domain or part of it
	 * \param v_size grid size on this physical domain
	 *
	 * \return An std::vector representing the local grid on each dimension
	 *
	 */
	std::vector<size_t> getGridSize(SpaceBox<dim,typename Decomposition::domain_type> & sp, Box<dim,typename Decomposition::domain_type> & domain, size_t (& v_size)[dim])
	{
		std::vector<size_t> tmp;
		for (size_t d = 0 ; d < dim ; d++)
		{
			//! Get the grid size compared to the domain space and its resolution
			typename Decomposition::domain_type dim_sz = (sp.getHigh(d) - sp.getLow(d)) / ((domain.getHigh(d) - domain.getLow(d)) / v_size[d]) + 0.5;

			// push the size of the local grid
			tmp.push_back(dim_sz);
		}
		return tmp;
	}

	/*! \brief Get the grid size
	 *
	 * Get the grid size, given a spaceBox
	 * it give the size on all directions of the local grid
	 *
	 * \param sp SpaceBox enclosing the local grid
	 * \param sz array to fill with the local grid size on each dimension
	 *
	 */
	void getGridSize(SpaceBox<dim,size_t> & sp, size_t (& v_size)[dim])
	{
		for (size_t d = 0 ; d < dim ; d++)
		{
			// push the size of the local grid
			v_size[d] = sp.getHigh(d) - sp.getLow(d);
		}
	}

//	bool link_ig_eg_init = false;

	/*! \brief Link the internal ghost boxes with the internal ghost boxes
	 *
	 * Each internal ghost box is linked with the external ghost box of the neighbor
	 * processor, in this function, the processors send the external ghost boxes to
	 * the near processors and link the internal ghost boxes to the received external
	 *
	 */
/*	void link_ig_eg()
	{
		if (link_ig_eg == true)	return;

		openfpm::vector< openfpm::vector< ::Box<dim,T>> > e_box_link(eg_box.size());
		openfpm::vector<size_t> prc_b;

		// Create a vector with external ghost boxes to send for each processor
		for (size_t i = 0; i < eg_box.size() ; i++)
		{
			for (size_t j = 0 ; j < eg_box.get(i).bid.size() ;j++)
			{
				e_box_link.add(eg_box.get(i).bid.get(j));
			}

			prc_b.add(eg_box.get(i).prc);
		}

		// Exchange the information

		v_cl.sendrecvMultipleMessagesNBX(prc_b,e_box_link,msg_alloc_external_box,this);

		// create a vector of boxes from the received messages
		//! None eg_box.size() == ig_box.size() == dec.getNNProcessors()
		for (size_t i = 0; i < dec.getNNProcessors() ; i++)
		{
			size_t n_ele = recv_sz.get(i) / sizeof(::Box<dim,T>);

			// Pointer of the received positions for each near processor
			void * ptr_boxes = recv_mem_gg.get(i).getPointer();

			PtrMemory * ptr1 = new PtrMemory(ptr_boxes,n_ele * sizeof(point));

			// received external ghost boxes in vector representation
			openfpm::vector< ::Box<dim,T>,openfpm::device_cpu< ::Box<dim,T>> ,PtrMemory,openfpm::grow_policy_identity> r_eg_box;

			// for each received external boxes
			for (size_t j = 0 ; j < r_eg_box.size() ; j++)
			{
				// get the middle point
				Point<dim,T> CM = r_eg_box.get(j).middle();

				// internal ghost with maximum volume
				long int max_vol_id = -1;
				T max_vol = 0;

				// Get the internal ghost boxes that fall into this point
				auto b_it = dec.getInternalIDBoxes();

				// Here we intersect each received external box with all
				// our internal ghost boxes, in theory only one box should match
				// the intersection, but we take the one with maximum, volume intersection
				while (b_it.isNext())
				{
					size_t b_id = b_it.get();

					// Get the internal ghost box
					const Box<dim,T> & b = dec.getIGhostBox(b_id);

					// out intersection
					Box<dim,T> b_out;
					// intersect
					bool intersect = b.Intersect(b,b_out);

					// if intersect
					if (intersect == true && b_out.getVolume() > max_vol)
					{
						max_vol = b_out.getVolume();
						max_vol_id = b_id;
					}

					++b_it;
				}

				// Link

				ig_box.get(i).bid.get(max_vol_id).r_id = j;
			}
		}

		link_ig_eg_init = true;
	}*/

	// Receiving size
	openfpm::vector<size_t> recv_sz;

	// Receiving buffer for particles ghost get
	openfpm::vector<HeapMemory> recv_mem_gg;

	/*! \brief Call-back to allocate buffer to receive incoming objects (external ghost boxes)
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
	static void * msg_alloc_external_box(size_t msg_i ,size_t total_msg, size_t total_p, size_t i, size_t ri, void * ptr)
	{
		grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> * g = static_cast<grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> *>(ptr);

		g->recv_sz.resize(g->dec.getNNProcessors());
		g->recv_mem_gg.resize(g->dec.getNNProcessors());

		// Get the local processor id
		size_t lc_id = g->dec.ProctoID(i);

		// resize the receive buffer
		g->recv_mem_gg.get(lc_id).resize(msg_i);
		g->recv_sz.get(lc_id) = msg_i;

		return g->recv_mem_gg.get(lc_id).getPointer();
	}

	/*! \brief Create per-processor internal ghost box list in grid units
	 *
	 * It also create g_id_to_external_ghost_box
	 *
	 */
	void create_ig_box()
	{
		// Get the grid info
		auto g = cd_sm.getGrid();

		if (init_i_g_box == true)	return;

		// Get the number of near processors
		for (size_t i = 0 ; i < dec.getNNProcessors() ; i++)
		{
			ig_box.add();
			auto&& pib = ig_box.last();

			pib.prc = dec.IDtoProc(i);
			for (size_t j = 0 ; j < dec.getProcessorNIGhost(i) ; j++)
			{
				// Get the internal ghost boxes and transform into grid units
				::Box<dim,St> ib = dec.getProcessorIGhostBox(i,j);
				ib /= cd_sm.getCellBox().getP2();

				// save the box and the sub-domain id (it is calculated as the linearization of P1)
				// It is unique because it is ensured that boxes does not overlap
				::Box<dim,size_t> cvt = ib;

				i_box_id bid_t;
				bid_t.box = cvt;
				bid_t.g_id = g.LinId(bid_t.box.getKP1());
				bid_t.sub = dec.getProcessorIGhostSub(i,j);
				pib.bid.add(bid_t);

				// Add the element in the unordered map
				g_id_to_external_ghost_box[bid_t.g_id] = bid_t.sub;
			}
		}

		init_i_g_box = true;
	}

	/*! \brief Create per-processor internal ghost box list in grid units
	 *
	 */
	void create_eg_box()
	{
		// Get the grid info
		auto g = cd_sm.getGrid();

		if (init_e_g_box == true)	return;

		// Get the number of near processors
		for (size_t i = 0 ; i < dec.getNNProcessors() ; i++)
		{
			eg_box.add();
			auto&& pib = eg_box.last();

			pib.prc = dec.IDtoProc(i);
			for (size_t j = 0 ; j < dec.getProcessorNEGhost(i) ; j++)
			{
				// Get the internal ghost boxes and transform into grid units
				::Box<dim,St> ib = dec.getProcessorEGhostBox(i,j);
				ib /= cd_sm.getCellBox().getP2();

				// save the box and the unique external ghost box id (linearization of P1)
				// It is (locally) unique because it is ensured that external ghost boxes does not overlap
				::Box<dim,size_t> cvt = ib;

				// sub domain id at which belong the external ghost box
				size_t sub_id = dec.getProcessorEGhostSub(i,j);

				e_box_id bid_t;
				bid_t.sub = sub_id;
				bid_t.g_e_box = cvt;
				bid_t.l_e_box = cvt;
				// Translate in local coordinate
				bid_t.l_e_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

				pib.bid.add(bid_t);
			}
		}

		init_e_g_box = true;
	}

public:

	//! constructor
	grid_dist_id(Vcluster v_cl, Decomposition & dec, const size_t (& g_sz)[dim], const Box<dim,St> & domain, const Ghost<dim,T> & ghost)
	:domain(domain),ghost(ghost),loc_grid(NULL),cd_sm(domain,g_sz,0),v_cl(v_cl),dec(dec)
	{
		// fill the global size of the grid
		for (int i = 0 ; i < dim ; i++)	{this->g_sz[i] = g_sz[i];}

		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * SUB_UNIT_FACTOR;

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[dim];
		for (int i = 0 ; i < dim ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/dim));}

		// Create the sub-domains
		dec.setParameters(div);

		// Create local grid
		Create();

		// Calculate ghost boxes
		dec.calculateGhostBoxes(ghost);
	}

	/*! \brief Constrcuctor
	 *
	 * \param g_sz array with the grid size on each dimension
	 * \param domain
	 *
	 */
	grid_dist_id(const size_t (& g_sz)[dim],const Box<dim,St> & domain, const Ghost<dim,St> & g)
	:domain(domain),ghost(g),dec(Decomposition(*global_v_cluster)),cd_sm(domain,g_sz,0),v_cl(*global_v_cluster)
	{
		// fill the global size of the grid
		for (size_t i = 0 ; i < dim ; i++)	{this->g_sz[i] = g_sz[i];}

		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * SUB_UNIT_FACTOR;

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/dim));}

		// Create the sub-domains
		dec.setParameters(div,domain);

		// Create local grid
		Create();

		// Calculate ghost boxes
		dec.calculateGhostBoxes(ghost);
	}

	/*! \brief Get the object that store the decomposition information
	 *
	 * \return the decomposition object
	 *
	 */

	Decomposition & getDecomposition()
	{
		return dec;
	}

	/*! \brief Create the grid on memory
	 *
	 */

	void Create()
	{
		// Box used for rounding error
		Box<dim,St> rnd_box;
		for (size_t i = 0 ; i < dim ; i++)	{rnd_box.setHigh(i,0.5); rnd_box.setLow(i,0.5);}
		// Box used for rounding in case of ghost
		Box<dim,St> g_rnd_box;
		for (size_t i = 0 ; i < dim ; i++)	{g_rnd_box.setHigh(i,0.5); g_rnd_box.setLow(i,-0.5);}

		// ! Create an hyper-cube approximation.
		// ! In order to work on grid_dist the decomposition
		// ! has to be a set of hyper-cube
		dec.hyperCube();

		// Get the number of local grid needed
		size_t n_grid = dec.getNLocalHyperCube();

		// create local grids for each hyper-cube
		loc_grid = v_cl.allocate<device_grid>(n_grid);

		// Size of the grid on each dimension
		size_t l_res[dim];

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{
			gdb_ext.add();

			// Get the local hyper-cube
			SpaceBox<dim,St> sp = dec.getLocalHyperCube(i);

			// Convert sp into grid units
			sp /= cd_sm.getCellBox().getP2();

			// enlarge by 0.5 for rounding
			sp.enlarge(rnd_box);

			// Convert from SpaceBox<dim,float> to SpaceBox<dim,long int>
			SpaceBox<dim,long int> sp_t = sp;

			// convert the ghost from space coordinate to grid units
			Ghost<dim,St> g_int = ghost;
			g_int /= cd_sm.getCellBox().getP2();

			// enlarge by 0.5 for rounding
			g_int.enlarge(g_rnd_box);

			// convert from Ghost<dim,St> to Ghost<dim,long int>
			Ghost<dim,long int> g_int_t = g_int;

			//! Save the origin of the local grid
			gdb_ext.last().origin = sp_t.getP1();

			// Center the local grid to zero
			sp_t -= sp_t.getP1();

			// save information about the local grid: domain box seen inside the domain + ghost box (see GDBoxes for a visual meaning)
			// and where the GDBox start, or the origin of the local grid (+ghost) in global coordinate
			gdb_ext.last().Dbox = sp_t;
			gdb_ext.last().Dbox -= g_int_t.getP1();
			// needed because the last key coordinate is size - 1 on each direction
			gdb_ext.last().Dbox.shrinkP2(1);
			// The origin is the Domain box + ghost, so shift
			gdb_ext.last().origin += g_int_t.getP1();

			// Enlarge sp with the Ghost size
			sp_t.enlarge_fix_P1(g_int_t);

			// Get the size of the local grid
			for (size_t i = 0 ; i < dim ; i++) {l_res[i] = sp_t.getHigh(i);}

			// Set the dimensions of the local grid
			loc_grid.get(i).template resize<Memory>(l_res);
		}
	}

	/*! \brief It return an iterator of the bulk part of the grid with a specified margin
	 *
	 * For margin we mean that every point is at least m points far from the border
	 *
	 * \param m margin
	 *
	 * \return An iterator to a grid with specified margins
	 *
	 */
	grid_dist_iterator<dim,device_grid> getDomainIterator()
	{
		grid_dist_iterator<dim,device_grid> it(loc_grid,gdb_ext);

		return it;
	}

	//! Destructor
	~grid_dist_id()
	{
	}

	/*! \brief Get the Virtual Cluster machine
	 *
	 * \return the Virtual cluster machine
	 *
	 */

	Vcluster & getVC()
	{
		return v_cl;
	}

	/*! \brief Get the reference of the selected element
	 *
	 *
	 * \param p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 */
	template <unsigned int p>inline auto get(grid_dist_key_dx<dim> & v1) -> typename std::add_lvalue_reference<decltype(loc_grid.get(v1.getSub()).template get<p>(v1.getKey()))>::type
	{
		return loc_grid.get(v1.getSub()).template get<p>(v1.getKey());
	}

	/*! \brief it store a box, its unique id and the sub-domain from where it come from
	 *
	 */
	struct i_box_id
	{
		//! Box
		::Box<dim,size_t> box;

		//! id
		size_t g_id;

		//! sub
		size_t sub;
	};

	/*! \brief It store the information about the external ghost box
	 *
	 *
	 */
	struct e_box_id
	{
		//! Box defining the external ghost box in global coordinates
		::Box<dim,size_t> g_e_box;

		//! Box defining the external ghost box in local coordinates
		::Box<dim,size_t> l_e_box;

		//! sub_id in which sub-domain this box live
		size_t sub;
	};

	/*! \brief Per-processor Internal ghost box
	 *
	 */
	struct ip_box_grid
	{
		// ghost in grid units
		openfpm::vector<i_box_id> bid;

		//! processor id
		size_t prc;
	};

	/*! \brief Per-processor external ghost box
	 *
	 */
	struct ep_box_grid
	{
		// ghost in grid units
		openfpm::vector<e_box_id> bid;

		//! processor id
		size_t prc;
	};

	//! Memory for the ghost sending buffer
	Memory g_send_prp_mem;

	//! Memory for the ghost sending buffer
	Memory g_recv_prp_mem;

	//! Flag that indicate if the external ghost box has been initialized
	bool init_e_g_box = false;

	//! Flag that indicate if the internal ghost box has been initialized
	bool init_i_g_box = false;

	//! Internal ghost boxes in grid units
	openfpm::vector<ip_box_grid> ig_box;

	//! External ghost boxes in grid units
	openfpm::vector<ep_box_grid> eg_box;

	/*! \brief It synchronize getting the ghost part of the grid
	 *
	 * \tparam prp Properties to get (sequence of properties ids)
	 * \opt options (unused)
	 *
	 */
	template<int... prp> void ghost_get()
	{
		// Sending property object
		typedef object<typename object_creator<typename T::type,prp...>::type> prp_object;

		// Convert the ghost  internal boxes into grid unit boxes
		create_ig_box();

		// Convert the ghost external boxes into grid unit boxes
		create_eg_box();

		// total number of sending vector
		std::vector<size_t> pap_prp;

		// Create a packing request vector
		for ( size_t i = 0 ; i < ig_box.size() ; i++ )
		{
			// for each ghost box
			for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
			{
				// And linked sub-domain
				size_t sub_id = ig_box.get(i).bid.get(j).sub;
				// Internal ghost box
				Box<dim,size_t> g_ig_box = ig_box.get(i).bid.get(j).box;
				g_ig_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();

				// Pack a size_t for the internal ghost id
				Packer<size_t,HeapMemory>::packRequest(pap_prp);
				// Create a sub grid iterator spanning the internal ghost layer
				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_ig_box.getKP1(),g_ig_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template packRequest<prp...>(loc_grid.get(sub_id),sub_it,pap_prp);
			}
		}

		// resize the property buffer memory
		g_send_prp_mem.resize(ExtPreAlloc<Memory>::calculateMem(pap_prp));

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prAlloc_prp = *(new ExtPreAlloc<Memory>(pap_prp,g_send_prp_mem));
		prAlloc_prp.incRef();

		// Pack information
		Pack_stat sts;

		// Pack the information for each processor and send it
		for ( size_t i = 0 ; i < ig_box.size() ; i++ )
		{
			sts.mark();

			// for each ghost box
			for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
			{
				// And linked sub-domain
				size_t sub_id = ig_box.get(i).bid.get(j).sub;
				// Internal ghost box
				Box<dim,size_t> g_ig_box = ig_box.get(i).bid.get(j).box;
				g_ig_box -= gdb_ext.get(sub_id).origin.template convertPoint<size_t>();
				// Ghost box global id
				size_t g_id = ig_box.get(i).bid.get(j).g_id;

				// Pack a size_t for the internal ghost id
				Packer<size_t,HeapMemory>::pack(prAlloc_prp,g_id,sts);
				// Create a sub grid iterator spanning the internal ghost layer
				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_ig_box.getKP1(),g_ig_box.getKP2());
				// and pack the internal ghost grid
				Packer<device_grid,HeapMemory>::template pack<prp...>(prAlloc_prp,loc_grid.get(sub_id),sub_it,sts);
			}
			// send the request
			v_cl.send(ig_box.get(i).prc,0,sts.getMarkPointer(prAlloc_prp),sts.getMarkSize(prAlloc_prp));
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

				//
				prp_recv[prp_recv.size()-1] += g_eg_box.getVolumeKey() * sizeof(prp_object) + sizeof(size_t);
			}
		}

		//! Resize the receiving buffer
		g_recv_prp_mem.resize(ExtPreAlloc<Memory>::calculateMem(prp_recv));

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prRecv_prp = *(new ExtPreAlloc<Memory>(prp_recv,g_recv_prp_mem));
		prRecv_prp.incRef();

		// queue the receive

		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			v_cl.recv(eg_box.get(i).prc,0,prRecv_prp.getPointer(i),prp_recv[i]);
		}

		// wait to receive communication
		v_cl.execute();

		Pack_stat ps;

		// Unpack the object
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			// for each external ghost box
			for (size_t j = 0 ; j < eg_box.get(i).bid.size() ; j++)
			{
				// Unpack the ghost box global-id

				size_t g_id;
				Unpacker<size_t,HeapMemory>::unpack(prRecv_prp,g_id,ps);

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
				Unpacker<device_grid,HeapMemory>::template unpack<prp...>(prRecv_prp,sub2,loc_grid.get(sub_id),ps);
			}
		}
	}

	/*! \brief Write the grid_dist_id information as VTK file
	 *
	 * The function generate several files
	 *
	 * 1)
	 *
	 * where X is the processor number
	 *
	 * \param output directory where to write the files
	 *
	 */
	bool write(std::string output) const
	{
		VTKWriter<openfpm::vector<::Box<dim,size_t>>,VECTOR_BOX> vtk_box1;

		openfpm::vector< openfpm::vector< ::Box<dim,size_t> > > boxes;

		//! Carefully we have to ensure that boxes does not reallocate inside the for loop
		boxes.reserve(ig_box.size());

		//! Write internal ghost in grid units (Color encoded)
		for (size_t p = 0 ; p < ig_box.size() ; p++)
		{
			boxes.add();

			// Create a vector of boxes
			for (size_t j = 0 ; j < ig_box.get(p).bid.size() ; j++)
			{
				boxes.last().add(ig_box.get(p).bid.get(j).box);
			}

			vtk_box1.add(boxes.last());
		}
		vtk_box1.write(output + std::string("internal_ghost_") + std::to_string(v_cl.getProcessUnitID()) + std::string(".vtk"));

		return true;
	}
};

/*! \brief This is a distributed grid
 *
 * Implementation of a distributed grid with id decomposition. A distributed grid is a grid distributed
 * across processors. The decomposition is performed on the id of the elements
 *
 * 1D specialization
 *
 * \param dim Dimensionality of the grid
 * \param T type of grid
 * \param Decomposition Class that decompose the grid for example CartDecomposition
 * \param Mem Is the allocator
 * \param device type of base structure is going to store the data
 *
 */

template<typename T, typename Decomposition,typename Memory , typename device_grid >
class grid_dist_id<1,T,Decomposition,Memory,device_grid>
{
	// Ghost
	Ghost<1,T> ghost;

	//! Local grids
	Vcluster_object_array<device_grid> loc_grid;

	//! Size of the grid on each dimension
	size_t g_sz[1];

	//! Communicator class
	Vcluster & v_cl;

	//! Extension of each grid: Domain and ghost + domain
	openfpm::vector<GBoxes<device_grid::dims>> gdb_ext;

	/*! \brief Get the grid size
	 *
	 * Get the grid size, given a domain, the resolution on it and another spaceBox
	 * it give the size on all directions of the local grid
	 *
	 * \param sp SpaceBox enclosing the local grid
	 * \param domain Space box enclosing the physical domain or part of it
	 * \param v_size grid size on this physical domain
	 *
	 * \return An std::vector representing the local grid on each dimension
	 *
	 */
	std::vector<size_t> getGridSize(SpaceBox<1,typename Decomposition::domain_type> & sp, Box<1,typename Decomposition::domain_type> & domain, size_t (& v_size)[1])
	{
		std::vector<size_t> tmp;
		for (size_t d = 0 ; d < 1 ; d++)
		{
			// push the size of the local grid
			tmp.push_back(g_sz[0]);
		}
		return tmp;
	}

	/*! \brief Get the grid size
	 *
	 * Get the grid size, given a spaceBox
	 * it give the size on all directions of the local grid
	 *
	 * \param sp SpaceBox enclosing the local grid
	 * \param sz array to fill with the local grid size on each dimension
	 *
	 */
	void getGridSize(SpaceBox<1,size_t> & sp, size_t (& v_size)[1])
	{
		for (size_t d = 0 ; d < 1 ; d++)
		{
			// push the size of the local grid
			v_size[d] = sp.getHigh(d) - sp.getLow(d);
		}
	}

public:

	//! constructor
	grid_dist_id(Vcluster v_cl, Decomposition & dec, size_t (& g_sz)[1], Box<1,T> & ghost)
	:ghost(ghost),loc_grid(NULL),v_cl(v_cl)
	{
		// All this code is completely untested to assume broken
		std::cout << "Error: " << __FILE__ << ":" << __LINE__ << " this structure is untested to assume broken ";

		// fill the global size of the grid
		for (int i = 0 ; i < 1 ; i++)	{this->g_sz[i] = g_sz[i];}

		// Create local memory
		Create();
	}

	//! constructor
	grid_dist_id(size_t (& g_sz)[1])
	:v_cl(*global_v_cluster),ghost(0)
	{
		// All this code is completely untested to assume broken
		std::cout << "Error: " << __FILE__ << ":" << __LINE__ << " this structure is untested to assume broken ";

		// fill the global size of the grid
		for (int i = 0 ; i < 1 ; i++)	{this->g_sz[i] = g_sz[i];}

		// Create local memory
		Create();
	}

	/*! \brief Create the grid on memory
	 *
	 */

	void Create()
	{
		size_t n_grid = 1;

		// create local grids for each hyper-cube
		loc_grid = v_cl.allocate<device_grid>(n_grid);

		// Size of the grid on each dimension
		size_t l_res[1];

		// Calculate the local grid size
		l_res[0] = g_sz[0] / v_cl.getProcessingUnits();

		// Distribute the remaining
		size_t residual = g_sz[0] % v_cl.getProcessingUnits();
		if (v_cl.getProcessUnitID() < residual)
			l_res[0]++;

		// Set the dimensions of the local grid
		loc_grid.get(0).template resize<Memory>(l_res);
	}

	/*! \brief It return an iterator of the bulk part of the grid with a specified margin
	 *
	 * For margin we mean that every point is at least m points far from the border
	 *
	 * \param m margin
	 *
	 * \return An iterator to a grid with specified margins
	 *
	 */
	grid_dist_iterator<1,device_grid> getDomainIterator()
	{
		grid_dist_iterator<1,device_grid> it(loc_grid,gdb_ext);

		return it;
	}

	//! Destructor
	~grid_dist_id()
	{
	}

	/*! \brief Get the Virtual Cluster machine
	 *
	 * \return the Virtual cluster machine
	 *
	 */

	Vcluster & getVC()
	{
		return v_cl;
	}

	/*! \brief Get the reference of the selected element
	 *
	 *
	 * \param p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 */
	template <unsigned int p>inline auto get(grid_dist_key_dx<1> & v1) -> typename std::add_lvalue_reference<decltype(loc_grid.get(v1.getSub()).template get<p>(v1.getKey()))>::type
	{
		return loc_grid.get(0).template get<p>(v1.getKey());
	}
};

#endif
