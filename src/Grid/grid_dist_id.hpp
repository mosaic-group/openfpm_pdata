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
 * Implementation of a distributed grid the decomposition is geometrical, grid
 * is splitted across several processor
 *
 * \param dim Dimensionality of the grid
 * \param St Type of space where the grid is living
 * \param T object the grid is storing
 * \param Decomposition Class that decompose the grid for example CartDecomposition
 * \param Mem Is the allocator
 * \param device type of base structure is going to store the data
 *
 * ### Create a distributed grid and access it
 * \snippet grid_dist_id_unit_test.hpp Create and access a distributed grid
 * ### Synchronize the ghosts and check the information
 * \snippet grid_dist_id_unit_test.hpp Synchronize the ghost and check the information
 * ### Create and access a distributed grid for complex structures
 * \snippet grid_dist_id_unit_test.hpp Create and access a distributed grid complex
 * ### Synchronize a distributed grid for complex structures
 * \snippet grid_dist_id_unit_test.hpp Synchronized distributed grid complex
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

	/*! \brief Create per-processor internal ghost boxes list in grid units and g_id_to_external_ghost_box
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
				::Box<dim,St> ib_dom = dec.getProcessorIGhostBox(i,j);
				::Box<dim,size_t> ib = cd_sm.convertDomainSpaceIntoGridUnits(ib_dom);

				// Check if ib is valid if not it mean that the internal ghost does not contain information so skip it
				if (ib.isValid() == false)
					continue;

				// save the box and the sub-domain id (it is calculated as the linearization of P1)
				::Box<dim,size_t> cvt = ib;

				i_box_id bid_t;
				bid_t.box = cvt;
				bid_t.g_id = dec.getProcessorIGhostId(i,j);
				bid_t.sub = dec.getProcessorIGhostSub(i,j);
				pib.bid.add(bid_t);
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
				// Get the external ghost boxes and transform into grid units
				::Box<dim,St> ib_dom = dec.getProcessorEGhostBox(i,j);
				::Box<dim,size_t> ib = cd_sm.convertDomainSpaceIntoGridUnits(ib_dom);

				// Check if ib is valid if not it mean that the internal ghost does not contain information so skip it
				if (ib.isValid() == false)
					continue;

				// save the box and the unique external ghost box id (linearization of P1)
				// It is (locally) unique because it is ensured that external ghost boxes does not overlap
				// Carefull it is not unique from the internal ghost box

				// sub domain id at which belong the external ghost box
				size_t sub_id = dec.getProcessorEGhostSub(i,j);

				e_box_id bid_t;
				bid_t.sub = sub_id;
				bid_t.g_e_box = ib;
				bid_t.l_e_box = ib;
				// Translate in local coordinate
				Box<dim,long int> tb = ib;
				tb -= gdb_ext.get(sub_id).origin;
				bid_t.l_e_box = tb;

				pib.bid.add(bid_t);

				// Add the map between the global ghost box id and id of the external box in the vector
				size_t g_id = dec.getProcessorEGhostId(i,j);
				g_id_to_external_ghost_box[g_id] = pib.bid.size()-1;
			}
		}

		init_e_g_box = true;
	}

	bool init_local_i_g_box = false;

	/*! \brief Create local internal ghost box in grid units
	 *
	 */
	void create_local_ig_box()
	{
		// Get the grid info
		auto g = cd_sm.getGrid();

		if (init_local_i_g_box == true)	return;

		// Get the number of near processors
		for (size_t i = 0 ; i < dec.getNLocalHyperCube() ; i++)
		{
			loc_ig_box.add();
			auto&& pib = loc_ig_box.last();

			for (size_t j = 0 ; j < dec.getLocalNIGhost(i) ; j++)
			{
				// Get the internal ghost boxes and transform into grid units
				::Box<dim,St> ib_dom = dec.getLocalIGhostBox(i,j);
				::Box<dim,size_t> ib = cd_sm.convertDomainSpaceIntoGridUnits(ib_dom);

				// Check if ib is valid if not it mean that the internal ghost does not contain information so skip it
				if (ib.isValid() == false)
					continue;

				pib.bid.add();
				pib.bid.last().box = ib;
				pib.bid.last().sub = dec.getLocalIGhostSub(i,j);
				pib.bid.last().k = dec.getLocalIGhostE(i,j);
			}
		}

		init_local_i_g_box = true;
	}

	bool init_local_e_g_box = false;

	/*! \brief Create per-processor external ghost boxes list in grid units
	 *
	 */
	void create_local_eg_box()
	{
		// Get the grid info
		auto g = cd_sm.getGrid();

		if (init_local_e_g_box == true)	return;

		// Get the number of near processors
		for (size_t i = 0 ; i < dec.getNLocalHyperCube() ; i++)
		{
			loc_eg_box.add();
			auto&& pib = loc_eg_box.last();

			for (size_t j = 0 ; j < dec.getLocalNEGhost(i) ; j++)
			{
				// Get the internal ghost boxes and transform into grid units
				::Box<dim,St> ib_dom = dec.getLocalEGhostBox(i,j);
				::Box<dim,size_t> ib = cd_sm.convertDomainSpaceIntoGridUnits(ib_dom);

				// Warning even if the ib is not a valid in grid unit we are forced to keep it
				// otherwise the value returned from dec.getLocalEGhostSub(i,j) will point to an
				// invalid or wrong box

				pib.bid.add();
				pib.bid.last().box = ib;
				pib.bid.last().sub = dec.getLocalEGhostSub(i,j);
			}
		}

		init_local_e_g_box = true;
	}

	/*! \brief Sync the local ghost part
	 *
	 * \tparam prp... properties to sync
	 *
	 */
	template<int... prp> void ghost_get_local()
	{
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

				// local external ghost box connected
				size_t k = loc_ig_box.get(i).bid.get(j).k;

				Box<dim,size_t> bx_dst = loc_eg_box.get(sub_id_dst).bid.get(k).box;

				// convert into local
				bx_dst -= gdb_ext.get(sub_id_dst).origin;

				// create 2 sub grid iterator
				grid_key_dx_iterator_sub<dim> sub_src(loc_grid.get(i).getGrid(),bx_src.getKP1(),bx_src.getKP2());
				grid_key_dx_iterator_sub<dim> sub_dst(loc_grid.get(sub_id_dst).getGrid(),bx_dst.getKP1(),bx_dst.getKP2());

#ifdef DEBUG

				if (loc_eg_box.get(sub_id_dst).bid.get(k).sub != i)
					std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination are not correctly linked" << "\n";

				if (sub_src.getVolume() != sub_dst.getVolume())
					std::cerr << "Error " << __FILE__ << ":" << __LINE__ << " source and destination does not match in size" << "\n";

#endif

				const auto & gs = loc_grid.get(i);
				auto & gd = loc_grid.get(sub_id_dst);

				while (sub_src.isNext())
				{
					// Option 1
					gd.set(sub_dst.get(),gs,sub_src.get());

					// Option 2
//					gd.get_o(sub_dst.get()) = gs.get_o(sub_src.get());

					++sub_src;
					++sub_dst;
				}
			}
		}
	}

	/*! \brief Check the grid has a valid size
	 *
	 */
	inline void check_size(const size_t (& g_sz)[dim])
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (g_sz[i] < 2)
				std::cerr << "Error: " << __FILE__ << ":" << __LINE__ << " distrobuted grids with size smaller than 2 are not supported\n";
		}
	}


	void write_ie_boxes(std::string output)
	{
		// Write internal ghost box
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
	}

protected:


	/*! \brief Get the point where it start the origin of the grid in the sub-domain i
	 *
	 * \return the point
	 *
	 */
	Point<dim,St> getOffset(size_t i)
	{
		return Point<dim,St>(gdb_ext.get(i).origin) * cd_sm.getCellBox().getP2();
	}

	/*! \brief Given a local sub-domain i with a local grid Domain + ghost return the part of the local grid that is domain
	 *
	 * \return the Box defining the domain in the local grid
	 *
	 */
	Box<dim,size_t> getDomain(size_t i)
	{
		return gdb_ext.get(i).Dbox;
	}

public:

	typedef device_grid loc_grid_type;

	/*! \brief Constrcuctor
	 *
	 * \param g_sz array with the grid size on each dimension
	 * \param domain domain where this grid live
	 * \param g Ghost
	 *
	 */
	grid_dist_id(const size_t (& g_sz)[dim],const Box<dim,St> & domain, const Ghost<dim,St> & g)
	:domain(domain),ghost(g),dec(Decomposition(*global_v_cluster)),v_cl(*global_v_cluster)
	{
		// check that the grid has valid size
		check_size(g_sz);

		// For a 5x5 grid you have 4x4 Cell
		size_t c_g[dim];
		for (size_t i = 0 ; i < dim ; i++)	{c_g[i] = (g_sz[i]-1 > 0)?(g_sz[i]-1):1;}

		// Initialize the cell decomposer
		cd_sm.setDimensions(domain,c_g,0);

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
		dec.setParameters(div,domain,ghost);

		// Create local grid
		Create();

		// Calculate ghost boxes
		dec.calculateGhostBoxes();
	}

	/*! \brief Get the object that store the information about the decomposition
	 *
	 * \return the decomposition object
	 *
	 */
	Decomposition & getDecomposition()
	{
		return dec;
	}

	/*! \brief Return the cell decomposer
	 *
	 * \return the cell decomposer
	 *
	 */
	const CellDecomposer_sm<dim,St> & getCellDecomposer()
	{
		return cd_sm;
	}

	/*! \brief Create the grids on memory
	 *
	 */
	void Create()
	{
		Box<dim,St> g_rnd_box;
		for (size_t i = 0 ; i < dim ; i++)	{g_rnd_box.setHigh(i,0.5); g_rnd_box.setLow(i,-0.5);}

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
			SpaceBox<dim,St> sp_g = dec.getSubDomainWithGhost(i);

			// Convert from SpaceBox<dim,St> to SpaceBox<dim,long int>
			SpaceBox<dim,long int> sp_t = cd_sm.convertDomainSpaceIntoGridUnits(sp);
			SpaceBox<dim,long int> sp_tg = cd_sm.convertDomainSpaceIntoGridUnits(sp_g);

			//! Save the origin of the sub-domain of the local grid
			gdb_ext.last().origin = sp_tg.getP1();

			// save information about the local grid: domain box seen inside the domain + ghost box (see GDBoxes for a visual meaning)
			// and where the GDBox start, or the origin of the local grid (+ghost) in global coordinate
			gdb_ext.last().Dbox = sp_t;
			gdb_ext.last().Dbox -= sp_tg.getP1();

			// center to zero
			sp_tg -= sp_tg.getP1();

			// Get the size of the local grid
			for (size_t i = 0 ; i < dim ; i++) {l_res[i] = (sp_tg.getHigh(i) >= 0)?(sp_tg.getHigh(i)+1):0;}

			// Set the dimensions of the local grid
			loc_grid.get(i).resize(l_res);
		}
	}

	/*! \brief Check that the global grid key is inside the grid domain
	 *
	 * \return true if is inside
	 *
	 */
	bool isInside(const grid_key_dx<dim> & gk) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (gk.get(i) < 0 || gk.get(i) >= (long int)g_sz[i])
				return false;
		}

		return true;
	}

	/*! \brief It return an iterator that span the full grid domain (each processor span its local domain)
	 *
	 * \return the iterator
	 *
	 */
	grid_dist_iterator<dim,device_grid,FREE> getDomainIterator()
	{
		grid_dist_iterator<dim,device_grid,FREE> it(loc_grid,gdb_ext);

		return it;
	}

	/*! \brief It return an iterator that span the grid domain + ghost part
	 *
	 *
	 */
	grid_dist_iterator<dim,device_grid,FIXED> getDomainGhostIterator()
	{
		grid_dist_iterator<dim,device_grid,FIXED> it(loc_grid,gdb_ext);

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

	/*! \brief it store an internal ghost box, the linked external ghost box and the sub-domain from where
	 *  it come from as internal ghost box
	 *
	 */
	struct i_lbox_id
	{
		//! Box
		::Box<dim,size_t> box;

		//! sub-domain id
		size_t sub;

		//! external box
		size_t k;
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

	/*! \brief It store the information about the external ghost box
	 *
	 *
	 */
	struct e_lbox_id
	{
		//! Box defining the external ghost box in local coordinates
		::Box<dim,size_t> box;

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

	/*! \brief local Internal ghost box
	 *
	 */
	struct i_lbox_grid
	{
		// ghost in grid units
		openfpm::vector<i_lbox_id> bid;
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

	/*! \brief Per-processor external ghost box
	 *
	 */
	struct e_lbox_grid
	{
		// ghost in grid units
		openfpm::vector<e_lbox_id> bid;
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

	//! Local internal ghost boxes in grid units
	openfpm::vector<i_lbox_grid> loc_ig_box;

	//! Local external ghost boxes in grid units
	openfpm::vector<e_lbox_grid> loc_eg_box;

	/*! \brief It synchronize the ghost parts
	 *
	 * \tparam prp... Properties to synchronize
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

		// Convert the local ghost internal boxes into grid unit boxes
		create_local_ig_box();

		// Convert the local external ghost boxes into grid unit boxes
		create_local_eg_box();

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
				prp_recv[prp_recv.size()-1] += g_eg_box.getVolumeKey() * sizeof(prp_object) + sizeof(size_t);
			}
		}

		//! Resize the receiving buffer
		g_recv_prp_mem.resize(ExtPreAlloc<Memory>::calculateMem(prp_recv));

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> & prRecv_prp = *(new ExtPreAlloc<Memory>(prp_recv,g_recv_prp_mem));
		prRecv_prp.incRef();

		// queue the receives
		for ( size_t i = 0 ; i < eg_box.size() ; i++ )
		{
			prRecv_prp.allocate(prp_recv[i]);
			v_cl.recv(eg_box.get(i).prc,0,prRecv_prp.getPointer(),prp_recv[i]);
		}

		// Before wait for the communication to complete we sync the local ghost
		// in order to overlap with communication

		ghost_get_local<prp...>();

		// wait to receive communication
		v_cl.execute();

		Unpack_stat ps;

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

	/*! \brief Get the spacing on each dimension
	 *
	 * \param get the spacing
	 *
	 */
	Point<dim,St> getSpacing()
	{
		return cd_sm.getCellBox().getP2();
	}

	/*! \brief Convert a g_dist_key_dx into a global key
	 *
	 * \see grid_dist_key_dx
	 * \see grid_dist_iterator
	 *
	 * \return the global position in the grid
	 *
	 */
	inline grid_key_dx<dim> getGKey(const grid_dist_key_dx<dim> & k)
	{
		// Get the sub-domain id
		size_t sub_id = k.getSub();

		grid_key_dx<dim> k_glob = k.getKey();

		// shift
		k_glob = k_glob + gdb_ext.get(sub_id).origin;

		return k_glob;
	}

	/*! \brief Write the distributed grid information
	 *
	 * * grid_X.vtk Output each local grids for each local processor X
	 * * internal_ghost_X.vtk Internal ghost boxes in grid units for the local processor X
	 *
	 * \param output directory where to put the files
	 *
	 */
	bool write(std::string output)
	{
		// Create a writer and write
		VTKWriter<boost::mpl::pair<device_grid,float>,VECTOR_GRIDS> vtk_g;
		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{
			Point<dim,St> offset = getOffset(i);
			vtk_g.add(loc_grid.get(i),offset,cd_sm.getCellBox().getP2(),gdb_ext.get(i).Dbox);
		}
		vtk_g.write(output + "/grid_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk");

		write_ie_boxes(output);

		return true;
	}

	/*! \brief Get the i sub-domain grid
	 *
	 * \param i sub-domain
	 *
	 * \return local grid
	 *
	 */
	device_grid & get_loc_grid(size_t i)
	{
		return loc_grid.get(i);
	}

	/*! \brief Return the number of local grid
	 *
	 * \return the number of local grid
	 *
	 */
	size_t getN_loc_grid()
	{
		return loc_grid.size();
	}
};



#endif
