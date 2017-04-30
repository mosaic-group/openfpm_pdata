#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include <unordered_map>
#include "Grid/map_grid.hpp"
#include "VCluster/VCluster.hpp"
#include "Space/SpaceBox.hpp"
#include "util/mathutil.hpp"
#include "Iterators/grid_dist_id_iterator_dec.hpp"
#include "Iterators/grid_dist_id_iterator.hpp"
#include "Iterators/grid_dist_id_iterator_sub.hpp"
#include "grid_dist_key.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "Packer_Unpacker/Packer.hpp"
#include "Packer_Unpacker/Unpacker.hpp"
#include "Decomposition/CartDecomposition.hpp"
#include "data_type/aggregate.hpp"

//! Internal ghost box sent to construct external ghost box into the other processors
template<unsigned int dim>
struct Box_fix
{
	//! Box in global unit
	Box<dim,size_t> bx;
	//! In which sector live the box
	comb<dim> cmb;
	//! Global id of the internal ghost box
	size_t g_id;
	//! from which sub-domain this internal ghost box is generated (or with which sub-domain is overlapping)
	size_t r_sub;
};

#define GRID_SUB_UNIT_FACTOR 64

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
 * \snippet grid_dist_id_unit_test.cpp Create and access a distributed grid
 * ### Synchronize the ghosts and check the information
 * \snippet grid_dist_id_unit_test.cpp Synchronize the ghost and check the information
 * ### Create and access a distributed grid for complex structures
 * \snippet grid_dist_id_unit_test.cpp Create and access a distributed grid complex
 * ### Synchronize a distributed grid for complex structures
 * \snippet grid_dist_id_unit_test.cpp Synchronized distributed grid complex
 * ### Usage of a grid dist iterator sub
 * \snippet grid_dist_id_unit_test.cpp Usage of a sub_grid iterator
 * ### Construct two grid with the same decomposition
 * \snippet grid_dist_id_unit_test.cpp Construct two grid with the same decomposition
 *
 */
template<unsigned int dim, typename St, typename T, typename Decomposition = CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_id
{
	//! Domain
	Box<dim,St> domain;

	//! Ghost expansion
	Ghost<dim,St> ghost;

	//! Local grids
	openfpm::vector<device_grid> loc_grid;

	//! Space Decomposition
	Decomposition dec;

	//! Extension of each grid: Domain and ghost + domain
	openfpm::vector<GBoxes<device_grid::dims>> gdb_ext;

	//! Size of the grid on each dimension
	size_t g_sz[dim];

	//! Structure that divide the space into cells
	CellDecomposer_sm<dim,St,shift<dim,St>> cd_sm;

	//! Communicator class
	Vcluster & v_cl;

	//! It map a global ghost id (g_id) to the external ghost box information
	//! It is unique across all the near processor
	std::unordered_map<size_t,size_t> g_id_to_external_ghost_box;

	//! It map a global ghost id (g_id) to the internal ghost box information
	//! (is unique for processor), it is not unique across all the near processor
	openfpm::vector<std::unordered_map<size_t,size_t>> g_id_to_internal_ghost_box;

	//! Receiving size
	openfpm::vector<size_t> recv_sz;

	//! Receiving buffer for particles ghost get
	openfpm::vector<HeapMemory> recv_mem_gg;

	//! Grid informations object
	grid_sm<dim,T> ginfo;

	//! Grid informations object without type
	grid_sm<dim,void> ginfo_v;

	//! Indicate if the local internal ghost box has been initialized
	bool init_local_i_g_box = false;

	//! Indicate if the local external ghost box has been initialized
	bool init_local_e_g_box = false;

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

	/*! \brief flip box just convert and internal ghost box into an external ghost box
	 *
	 * \param box to convert
	 * \param cmb sector position of the box
	 *
	 * \return the converted box
	 *
	 */
	Box<dim,long int> flip_box(const Box<dim,long int> & box, const comb<dim> & cmb)
	{
		Box<dim,long int> flp;

		for (size_t i = 0 ; i < dim; i++)
		{
			if (cmb[i] == 0)
			{
				flp.setLow(i,box.getLow(i));
				flp.setHigh(i,box.getHigh(i));
			}
			else if (cmb[i] == 1)
			{
				flp.setLow(i,box.getLow(i) + ginfo.size(i));
				flp.setHigh(i,box.getHigh(i) + ginfo.size(i));
			}
			else if (cmb[i] == -1)
			{
				flp.setLow(i,box.getLow(i) - ginfo.size(i));
				flp.setHigh(i,box.getHigh(i) - ginfo.size(i));
			}
		}

		return flp;
	}

	/*! \brief Create per-processor internal ghost boxes list in grid units and g_id_to_external_ghost_box
	 *
	 */
	void create_ig_box()
	{
		if (init_i_g_box == true)	return;

		// Get the grid info
		auto g = cd_sm.getGrid();

		g_id_to_internal_ghost_box.resize(dec.getNNProcessors());

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
				::Box<dim,long int> ib = cd_sm.convertDomainSpaceIntoGridUnits(ib_dom,dec.periodicity());

				// Check if ib is valid if not it mean that the internal ghost does not contain information so skip it
				if (ib.isValid() == false)
					continue;

				// save the box and the sub-domain id (it is calculated as the linearization of P1)
				::Box<dim,size_t> cvt = ib;

				i_box_id<dim> bid_t;
				bid_t.box = cvt;
				bid_t.g_id = dec.getProcessorIGhostId(i,j);
				bid_t.sub = dec.getProcessorIGhostSub(i,j);
				bid_t.cmb = dec.getProcessorIGhostPos(i,j);
				bid_t.r_sub = dec.getProcessorIGhostSSub(i,j);
				pib.bid.add(bid_t);

				g_id_to_internal_ghost_box.get(i)[bid_t.g_id] = pib.bid.size()-1;
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

		// Here we collect all the calculated internal ghost box in the sector different from 0 that this processor has

		openfpm::vector<size_t> prc;
		openfpm::vector<size_t> prc_recv;
		openfpm::vector<size_t> sz_recv;
		openfpm::vector<openfpm::vector<Box_fix<dim>>> box_int_send(dec.getNNProcessors());
		openfpm::vector<openfpm::vector<Box_fix<dim>>> box_int_recv;

		for(size_t i = 0 ; i < dec.getNNProcessors() ; i++)
		{
			for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
			{
				box_int_send.get(i).add();
				box_int_send.get(i).last().bx = ig_box.get(i).bid.get(j).box;
				box_int_send.get(i).last().g_id = ig_box.get(i).bid.get(j).g_id;
				box_int_send.get(i).last().r_sub = ig_box.get(i).bid.get(j).r_sub;
				box_int_send.get(i).last().cmb = ig_box.get(i).bid.get(j).cmb;
			}
			prc.add(dec.IDtoProc(i));
		}

		v_cl.SSendRecv(box_int_send,box_int_recv,prc,prc_recv,sz_recv);

		eg_box.resize(dec.getNNProcessors());

		for (size_t i = 0 ; i < eg_box.size() ; i++)
			eg_box.get(i).prc = dec.IDtoProc(i);

		for (size_t i = 0 ; i < box_int_recv.size() ; i++)
		{
			size_t p_id = dec.ProctoID(prc_recv.get(i));
			auto&& pib = eg_box.get(p_id);
			pib.prc = prc_recv.get(i);

			// For each received internal ghost box
			for (size_t j = 0 ; j < box_int_recv.get(i).size() ; j++)
			{
				size_t send_list_id = box_int_recv.get(i).get(j).r_sub;

				// Get the list of the sent sub-domains
				// and recover the id of the sub-domain from
				// the sent list
				const openfpm::vector<size_t> & s_sub = dec.getSentSubdomains(p_id);
				size_t sub_id = s_sub.get(send_list_id);

				e_box_id<dim> bid_t;
				bid_t.sub = sub_id;
				bid_t.cmb = box_int_recv.get(i).get(j).cmb;
				bid_t.cmb.sign_flip();
				::Box<dim,long int> ib = flip_box(box_int_recv.get(i).get(j).bx,box_int_recv.get(i).get(j).cmb);
				bid_t.g_e_box = ib;
				bid_t.g_id = box_int_recv.get(i).get(j).g_id;
				// Translate in local coordinate
				Box<dim,long int> tb = ib;
				tb -= gdb_ext.get(sub_id).origin;
				bid_t.l_e_box = tb;

				pib.bid.add(bid_t);

				g_id_to_external_ghost_box[bid_t.g_id] = pib.bid.size()-1;
			}
		}

		init_e_g_box = true;
	}

	/*! \brief Create local internal ghost box in grid units
	 *
	 */
	void create_local_ig_box()
	{
		// Get the grid info
		auto g = cd_sm.getGrid();

		if (init_local_i_g_box == true)	return;

		// Get the number of sub-domains
		for (size_t i = 0 ; i < dec.getNSubDomain() ; i++)
		{
			loc_ig_box.add();
			auto&& pib = loc_ig_box.last();

			for (size_t j = 0 ; j < dec.getLocalNIGhost(i) ; j++)
			{
				// Get the internal ghost boxes and transform into grid units
				::Box<dim,St> ib_dom = dec.getLocalIGhostBox(i,j);
				::Box<dim,long int> ib = cd_sm.convertDomainSpaceIntoGridUnits(ib_dom,dec.periodicity());

				// Check if ib is valid if not it mean that the internal ghost does not contain information so skip it
				if (ib.isValid() == false)
					continue;

				pib.bid.add();
				pib.bid.last().box = ib;
				pib.bid.last().sub = dec.getLocalIGhostSub(i,j);
				pib.bid.last().k = dec.getLocalIGhostE(i,j);
				pib.bid.last().cmb = dec.getLocalIGhostPos(i,j);
			}
		}

		init_local_i_g_box = true;
	}

	/*! \brief Create per-processor external ghost boxes list in grid units
	 *
	 */
	void create_local_eg_box()
	{
		// Get the grid info
		auto g = cd_sm.getGrid();

		if (init_local_e_g_box == true)	return;

		loc_eg_box.resize(dec.getNSubDomain());

		// Get the number of sub-domain
		for (size_t i = 0 ; i < dec.getNSubDomain() ; i++)
		{
			for (size_t j = 0 ; j < loc_ig_box.get(i).bid.size() ; j++)
			{
				size_t k = loc_ig_box.get(i).bid.get(j).sub;
				auto & pib = loc_eg_box.get(k);

				size_t s = loc_ig_box.get(i).bid.get(j).k;
				pib.bid.resize(dec.getLocalNEGhost(k));

				pib.bid.get(s).box = flip_box(loc_ig_box.get(i).bid.get(j).box,loc_ig_box.get(i).bid.get(j).cmb);
				pib.bid.get(s).sub = dec.getLocalEGhostSub(k,s);
				pib.bid.get(s).cmb = loc_ig_box.get(i).bid.get(j).cmb;
				pib.bid.get(s).cmb.sign_flip();
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

				// local internal ghost box connected
				size_t k = loc_ig_box.get(i).bid.get(j).k;

				Box<dim,size_t> bx_dst = loc_eg_box.get(sub_id_dst).bid.get(k).box;

				// convert into local
				bx_dst -= gdb_ext.get(sub_id_dst).origin;

				// create 2 sub grid iterator

				if (bx_dst.isValid() == false)
					continue;

				grid_key_dx_iterator_sub<dim> sub_src(loc_grid.get(i).getGrid(),bx_src.getKP1(),bx_src.getKP2());
				grid_key_dx_iterator_sub<dim> sub_dst(loc_grid.get(sub_id_dst).getGrid(),bx_dst.getKP1(),bx_dst.getKP2());

#ifdef SE_CLASS1

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

					++sub_src;
					++sub_dst;
				}
			}
		}
	}

	/*! \brief Check the grid has a valid size
	 *
	 * \param g_sz size of the grid
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

	/*! \brief Create the grids on memory
	 *
	 */
	void Create()
	{
		// Get the number of local grid needed
		size_t n_grid = dec.getNSubDomain();

		// create gdb
		create_gdb_ext<dim,Decomposition>(gdb_ext,dec,cd_sm);

		// create local grids for each hyper-cube
		loc_grid.resize(n_grid);

		// Size of the grid on each dimension
		size_t l_res[dim];

		// Allocate the grids
		for (size_t i = 0 ; i < n_grid ; i++)
		{

			SpaceBox<dim,long int> sp_tg = gdb_ext.get(i).GDbox;

			// Get the size of the local grid
			// The boxes indicate the extension of the index the size
			// is this extension +1
			// for example a 1D box (interval) from 0 to 3 in one dimension have
			// the points 0,1,2,3 = so a total of 4 points
			for (size_t j = 0 ; j < dim ; j++)
				l_res[j] = (sp_tg.getHigh(j) >= 0)?(sp_tg.getHigh(j)+1):0;

			// Set the dimensions of the local grid
			loc_grid.get(i).resize(l_res);
		}
	}

	/*! \brief Default Copy constructor on this class make no sense and is unsafe, this definition disable it
	 *
	 * \param g grid to copy
	 *
	 */
	grid_dist_id(const grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> & g)
	:v_cl(g.v_cl)
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif
	}

    /*! \brief Initialize the Cell decomposer of the grid enforcing perfect overlap of the cells
	 *
	 * \param cd_old the CellDecomposer we are trying to mach
	 * \param ext extension of the domain
	 *
	 */
	inline void InitializeCellDecomposer(const CellDecomposer_sm<dim,St,shift<dim,St>> & cd_old, const Box<dim,size_t> & ext)
	{
		// Initialize the cell decomposer
		cd_sm.setDimensions(cd_old,ext);
	}

    /*! \brief Initialize the Cell decomposer of the grid
	 *
	 * \param g_sz Size of the grid
	 * \param bc boundary conditions
	 *
	 */
	inline void InitializeCellDecomposer(const size_t (& g_sz)[dim], const size_t (& bc)[dim])
	{
		// check that the grid has valid size
		check_size(g_sz);

		// get the size of the cell decomposer
		size_t c_g[dim];
		getCellDecomposerPar<dim>(c_g,g_sz,bc);

		// Initialize the cell decomposer
		cd_sm.setDimensions(domain,c_g,0);
	}

	/*! \brief Initialize the grid
	 *
	 * \param g_sz Global size of the grid
	 * \param bc boundary conditions
	 *
	 */
	inline void InitializeDecomposition(const size_t (& g_sz)[dim], const size_t (& bc)[dim])
	{
		// fill the global size of the grid
		for (size_t i = 0 ; i < dim ; i++)	{this->g_sz[i] = g_sz[i];}

		// Get the number of processor and calculate the number of sub-domain
		// for decomposition
		size_t n_proc = v_cl.getProcessingUnits();
		size_t n_sub = n_proc * GRID_SUB_UNIT_FACTOR;

		// Calculate the maximum number (before merging) of sub-domain on
		// each dimension
		size_t div[dim];
		for (size_t i = 0 ; i < dim ; i++)
		{div[i] = openfpm::math::round_big_2(pow(n_sub,1.0/dim));}

		// Create the sub-domains
		dec.setParameters(div,domain,bc,ghost);
		dec.decompose();
	}

	/*! \brief Initialize the grid
	 *
	 * \param g_sz Global size of the grid
	 *
	 */
	inline void InitializeStructures(const size_t (& g_sz)[dim])
	{
		// fill the global size of the grid
		for (size_t i = 0 ; i < dim ; i++)	{this->g_sz[i] = g_sz[i];}

		// Create local grid
		Create();
	}

protected:

	/*! \brief Get the point where it start the origin of the grid of the sub-domain i
	 *
	 * \param i sub-domain
	 *
	 * \return the point
	 *
	 */
	Point<dim,St> getOffset(size_t i)
	{
		return pmul(Point<dim,St>(gdb_ext.get(i).origin), cd_sm.getCellBox().getP2());
	}

	/*! \brief Given a local sub-domain i with a local grid Domain + ghost return the part of the local grid that is domain
	 *
	 * \param i sub-domain
	 *
	 * \return the Box defining the domain in the local grid
	 *
	 */
	Box<dim,size_t> getDomain(size_t i)
	{
		return gdb_ext.get(i).Dbox;
	}

	/*! \brief Convert a ghost from grid point units into continus space
	 *
	 * \param gd Ghost in continuous space
	 * \param cd_sm CellDecomposer of the grid
	 *
	 * \return the ghost in continuous unit
	 *
	 */
	static inline Ghost<dim,float> convert_ghost(const Ghost<dim,long int> & gd, const CellDecomposer_sm<dim,St,shift<dim,St>> & cd_sm)
	{
		Ghost<dim,float> gc;

		// get the grid spacing
		Box<dim,St> sp = cd_sm.getCellBox();

		// enlarge 0.001 of the spacing
		sp.magnify_fix_P1(1.1);

		// set the ghost
		for (size_t i = 0 ; i < dim ; i++)
		{
			gc.setLow(i,-sp.getHigh(i));
			gc.setHigh(i,sp.getHigh(i));
		}

		return gc;
	}

public:

	//! Which kind of grid the structure store
	typedef device_grid d_grid;

	//! Decomposition used
	typedef Decomposition decomposition;

	//! value_type
	typedef T value_type;

	//! Type of space
	typedef St stype;

	//! Type of Memory
	typedef Memory memory_type;

	//! Type of device grid
	typedef device_grid device_grid_type;

	//! Number of dimensions
	static const unsigned int dims = dim;

	/*! \brief Get the domain where the grid is defined
	 *
	 * \return the domain of the grid
	 *
	 */
	inline const Box<dim,St> getDomain() const
	{
		return domain;
	}

    /*! \brief Get the spacing of the grid in direction i
     *
     * \param i dimension
     *
     * \return the spacing
     *
     */
    inline St spacing(size_t i) const
    {
    	return cd_sm.getCellBox().getHigh(i);
    }

	/*! \brief Return the total number of points in the grid
	 *
	 * \return number of points
	 *
	 */
	size_t size() const
	{
		return ginfo_v.size();
	}

	/*! \brief Return the total number of points in the grid
	 *
	 * \param i direction
	 *
	 * \return number of points on direction i
	 *
	 */
	size_t size(size_t i) const
	{
		return ginfo_v.size(i);
	}

	/*! \brief This constructor is special, it construct an expanded grid that perfectly overlap with the previous
	 *
	 * The key-word here is "perfectly overlap". Using the default constructor you could create
	 * something similar, but because of rounding-off error it can happen that it is not perfectly overlapping
	 *
	 * \param g previous grid
	 * \param gh Ghost part in grid units
	 * \param ext extension of the grid (must be positive on every direction)
	 *
	 */
	template<typename H> grid_dist_id(const grid_dist_id<dim,St,H,typename Decomposition::base_type,Memory,grid_cpu<dim,H>> & g, const Ghost<dim,long int> & gh, Box<dim,size_t> ext)
	:dec(create_vcluster()),v_cl(create_vcluster())
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif

		size_t ext_dim[dim];
		for (size_t i = 0 ; i < dim ; i++) {ext_dim[i] = g.getGridInfoVoid().size(i) + ext.getKP1().get(i) + ext.getKP2().get(i);}

		// Set the grid info of the extended grid
		ginfo.setDimensions(ext_dim);
		ginfo_v.setDimensions(ext_dim);

		InitializeCellDecomposer(g.getCellDecomposer(),ext);

		ghost = convert_ghost(gh,cd_sm);

		// Extend the grid by the extension part and calculate the domain

		for (size_t i = 0 ; i < dim ; i++)
		{
			g_sz[i] = g.size(i) + ext.getLow(i) + ext.getHigh(i);

			if (g.getDecomposition().periodicity(i) == NON_PERIODIC)
			{
				this->domain.setLow(i,g.getDomain().getLow(i) - ext.getLow(i) * g.spacing(i) - g.spacing(i) / 2.0);
				this->domain.setHigh(i,g.getDomain().getHigh(i) + ext.getHigh(i) * g.spacing(i) + g.spacing(i) / 2.0);
			}
			else
			{
				this->domain.setLow(i,g.getDomain().getLow(i) - ext.getLow(i) * g.spacing(i));
				this->domain.setHigh(i,g.getDomain().getHigh(i) + ext.getHigh(i) * g.spacing(i));
			}
		}

		dec.setParameters(g.getDecomposition(),ghost,this->domain);

		InitializeStructures(g.getGridInfoVoid().getSize());
	}

    /*! It constructs a grid of a specified size, defined on a specified Box space, forcing to follow a specified decomposition and with a specified ghost size
     *
     * \param dec Decomposition
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param ghost Ghost part
     *
     */
    grid_dist_id(const Decomposition & dec, const size_t (& g_sz)[dim], const Ghost<dim,St> & ghost)
    :domain(dec.getDomain()),ghost(ghost),dec(dec),v_cl(create_vcluster()),ginfo(g_sz),ginfo_v(g_sz)
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif

		InitializeCellDecomposer(g_sz,dec.periodicity());
		InitializeStructures(g_sz);
	}

    /*! It constructs a grid of a specified size, defined on a specified Box space, forcing to follow a specified decomposition and with a specified ghost size
     *
     * \param dec Decomposition
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param ghost Ghost part
     *
     */
    grid_dist_id(Decomposition && dec, const size_t (& g_sz)[dim], const Ghost<dim,St> & ghost)
    :domain(dec.getDomain()),ghost(ghost),dec(dec),ginfo(g_sz),ginfo_v(g_sz),v_cl(create_vcluster())
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif

		InitializeCellDecomposer(g_sz,dec.periodicity());
		InitializeStructures(g_sz);
	}

    /*! It constructs a grid of a specified size, defined on a specified Box space, forcing to follow a specified decomposition, and having a specified ghost size
     *
     * \param dec Decomposition
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param g Ghost part (given in grid units)
     *
     * \warning In very rare case the ghost part can be one point bigger than the one specified
     *
     */
	grid_dist_id(const Decomposition & dec, const size_t (& g_sz)[dim], const Ghost<dim,long int> & g)
	:domain(dec.getDomain()),dec(create_vcluster()),v_cl(create_vcluster()),ginfo(g_sz),ginfo_v(g_sz)
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif

		InitializeCellDecomposer(g_sz,dec.periodicity());

		ghost = convert_ghost(g,cd_sm);
		this->dec = dec.duplicate(ghost);

		// Initialize structures
		InitializeStructures(g_sz);
	}

    /*! It construct a grid of a specified size, defined on a specified Box space, forcing to follow a specified decomposition, and having a specified ghost size
     *
     * \param dec Decomposition
     * \param g_sz grid size on each dimension
     * \param g Ghost part (given in grid units)
     *
     * \warning In very rare case the ghost part can be one point bigger than the one specified
     *
     */
	grid_dist_id(Decomposition && dec, const size_t (& g_sz)[dim], const Ghost<dim,long int> & g)
	:domain(dec.getDomain()),dec(dec),v_cl(create_vcluster()),ginfo(g_sz),ginfo_v(g_sz)
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif
		InitializeCellDecomposer(g_sz,dec.periodicity());

		ghost = convert_ghost(g,cd_sm);

		// Initialize structures
		InitializeStructures(g_sz);
	}

    /*! It construct a grid of a specified size, defined on a specified Box space, and having a specified ghost size
     *
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param g Ghost part (given in grid units)
     *
     * \warning In very rare case the ghost part can be one point bigger than the one specified
     *
     */
	grid_dist_id(const size_t (& g_sz)[dim],const Box<dim,St> & domain, const Ghost<dim,St> & g)
	:grid_dist_id(g_sz,domain,g,create_non_periodic<dim>())
	{
	}

    /*! It construct a grid of a specified size, defined on a specified Box space, having a specified ghost size and periodicity
     *
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param g Ghost part of the domain (given in grid units)
     *
     * \warning In very rare case the ghost part can be one point bigger than the one specified
     *
     */
	grid_dist_id(const size_t (& g_sz)[dim],const Box<dim,St> & domain, const Ghost<dim,long int> & g)
	:grid_dist_id(g_sz,domain,g,create_non_periodic<dim>())
	{
	}

    /*! It construct a grid of a specified size, defined on a specified Box space, having a specified ghost size, and specified periodicity
     *
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param g Ghost part (given in grid units)
     * \param p Boundary conditions
     *
     * \warning In very rare case the ghost part can be one point bigger than the one specified
     *
     */
	grid_dist_id(const size_t (& g_sz)[dim],const Box<dim,St> & domain, const Ghost<dim,St> & g, const periodicity<dim> & p)
	:domain(domain),ghost(g),dec(create_vcluster()),v_cl(create_vcluster()),ginfo(g_sz),ginfo_v(g_sz)
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif

		InitializeCellDecomposer(g_sz,p.bc);
		InitializeDecomposition(g_sz, p.bc);
		InitializeStructures(g_sz);
	}

    /*! It construct a grid of a specified size, defined on a specified Box space, having a specified ghost size and periodicity
     *
     * \param g_sz grid size on each dimension
     * \param domain Box that contain the grid
     * \param g Ghost part of the domain (given in grid units)
     * \param p periodicity
     *
     * \warning In very rare case the ghost part can be one point bigger than the one specified
     *
     */
	grid_dist_id(const size_t (& g_sz)[dim],const Box<dim,St> & domain, const Ghost<dim,long int> & g, const periodicity<dim> & p)
	:domain(domain),dec(create_vcluster()),v_cl(create_vcluster()),ginfo(g_sz),ginfo_v(g_sz)
	{
#ifdef SE_CLASS2
		check_new(this,8,GRID_DIST_EVENT,4);
#endif
		InitializeCellDecomposer(g_sz,p.bc);

		ghost = convert_ghost(g,cd_sm);

		InitializeDecomposition(g_sz,p.bc);
		// Initialize structures
		InitializeStructures(g_sz);
	}

	/*! \brief Get an object containing the grid informations
	 *
	 * \return an information object about this grid
	 *
	 */
	const grid_sm<dim,T> & getGridInfo() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return ginfo;
	}

	/*! \brief Get an object containing the grid informations without type
	 *
	 * \return an information object about this grid
	 *
	 */
	const grid_sm<dim,void> & getGridInfoVoid() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return ginfo_v;
	}

	/*! \brief Get the object that store the information about the decomposition
	 *
	 * \return the decomposition object
	 *
	 */
	Decomposition & getDecomposition()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return dec;
	}

	/*! \brief Get the object that store the information about the decomposition
	 *
	 * \return the decomposition object
	 *
	 */
	const Decomposition & getDecomposition() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return dec;
	}

	/*! \brief Return the cell decomposer
	 *
	 * \return the cell decomposer
	 *
	 */
	const CellDecomposer_sm<dim,St,shift<dim,St>> & getCellDecomposer() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return cd_sm;
	}

	/*! \brief Check that the global grid key is inside the grid domain
	 *
	 * \param gk point to check
	 *
	 * \return true if is inside
	 *
	 */
	bool isInside(const grid_key_dx<dim> & gk) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (gk.get(i) < 0 || gk.get(i) >= (long int)g_sz[i])
				return false;
		}

		return true;
	}

	/*! \brief Get the size of local domain grids
	 *
	 * \return The size of the local domain
	 *
	 */
	size_t getLocalDomainSize() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		size_t total = 0;

		for (size_t i = 0 ; i < gdb_ext.size() ; i++)
		{
			total += gdb_ext.get(i).Dbox.getVolumeKey();
		}

		return total;
	}

	/*! \brief It return the informations about the local grids
	 *
	 * \return The information about the local grids
	 *
	 */
	const openfpm::vector<GBoxes<device_grid::dims>> & getLocalGridsInfo()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gdb_ext;
	}

	/*! \brief It return an iterator that span the full grid domain (each processor span its local domain)
	 *
	 * \return the iterator
	 *
	 */
	grid_dist_iterator<dim,device_grid,FREE> getDomainIterator() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif

		grid_key_dx<dim> stop(ginfo_v.getSize());
		grid_key_dx<dim> one;
		one.one();
		stop = stop - one;

		grid_dist_iterator<dim,device_grid,FREE> it(loc_grid,gdb_ext,stop);

		return it;
	}

	/*! \brief It return an iterator that span the grid domain + ghost part
	 *
	 * \return the iterator
	 *
	 */
	grid_dist_iterator<dim,device_grid,FIXED> getDomainGhostIterator() const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		grid_dist_iterator<dim,device_grid,FIXED> it(loc_grid,gdb_ext);

		return it;
	}

	/*! \brief It return an iterator that span the grid domain only in the specified
	 * part
	 *
	 * The key spanned are the one inside the box spanned by the start point and the end
	 * point included
	 *
	 * \param start point
	 * \param stop point
	 *
	 * \return the sub-domain iterator
	 *
	 */
	grid_dist_iterator_sub<dim,device_grid> getSubDomainIterator(const grid_key_dx<dim> & start, const grid_key_dx<dim> & stop) const
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		grid_dist_iterator_sub<dim,device_grid> it(start,stop,loc_grid,gdb_ext);

		return it;
	}

	/*! \brief It return an iterator that span the grid domain only in the specified
	 * part
	 *
	 * The key spanned are the one inside the box spanned by the start point and the end
	 * point included
	 *
	 * \param start point
	 * \param stop point
	 *
	 * \return an iterator on the sub-part of the grid
	 *
	 */
	grid_dist_iterator_sub<dim,device_grid> getSubDomainIterator(const long int (& start)[dim], const long int (& stop)[dim]) const
	{
		grid_dist_iterator_sub<dim,device_grid> it(grid_key_dx<dim>(start),grid_key_dx<dim>(stop),loc_grid,gdb_ext);

		return it;
	}

	//! Destructor
	~grid_dist_id()
	{
#ifdef SE_CLASS2
		check_delete(this);
#endif
		dec.decRef();
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

	/*! \brief Indicate that this grid is not staggered
	 *
	 * \return false
	 *
	 */
	bool is_staggered()
	{
		return false;
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto get(const grid_dist_key_dx<dim> & v1) const -> typename std::add_lvalue_reference<decltype(loc_grid.get(v1.getSub()).template get<p>(v1.getKey()))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return loc_grid.get(v1.getSub()).template get<p>(v1.getKey());
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto get(const grid_dist_key_dx<dim> & v1) -> typename std::add_lvalue_reference<decltype(loc_grid.get(v1.getSub()).template get<p>(v1.getKey()))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return loc_grid.get(v1.getSub()).template get<p>(v1.getKey());
	}

	//! Memory for the ghost sending buffer
	Memory g_send_prp_mem;

	//! Memory for the ghost sending buffer
	Memory g_recv_prp_mem;

	//! Flag that indicate if the external ghost box has been initialized
	bool init_e_g_box = false;

	//! Flag that indicate if the internal ghost box has been initialized
	bool init_i_g_box = false;

	//! Flag that indicate if the internal and external ghost box has been fixed
	bool init_fix_ie_g_box = false;

	//! Internal ghost boxes in grid units
	openfpm::vector<ip_box_grid<dim>> ig_box;

	//! External ghost boxes in grid units
	openfpm::vector<ep_box_grid<dim>> eg_box;

	//! Local internal ghost boxes in grid units
	openfpm::vector<i_lbox_grid<dim>> loc_ig_box;

	//! Local external ghost boxes in grid units
	openfpm::vector<e_lbox_grid<dim>> loc_eg_box;

	/*! \brief It synchronize the ghost parts
	 *
	 * \tparam prp... Properties to synchronize
	 *
	 */
	template<int... prp> void ghost_get()
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif

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


		size_t req = 0;

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
		ExtPreAlloc<Memory> & prAlloc_prp = *(new ExtPreAlloc<Memory>(req,g_send_prp_mem));

		prAlloc_prp.incRef();

		// Pack information
		Pack_stat sts;

		// Pack the information for each processor and send it
		for ( size_t i = 0 ; i < ig_box.size() ; i++ )
		{

			sts.mark();
			void * pointer = prAlloc_prp.getPointerEnd();

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
				Packer<size_t,HeapMemory>::pack(prAlloc_prp,g_id,sts);
				// Create a sub grid iterator spanning the internal ghost layer
				grid_key_dx_iterator_sub<dim> sub_it(loc_grid.get(sub_id).getGrid(),g_ig_box.getKP1(),g_ig_box.getKP2());
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
		ExtPreAlloc<Memory> & prRecv_prp = *(new ExtPreAlloc<Memory>(tot_recv,g_recv_prp_mem));
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

	/*! \brief Copy the give grid into this grid
	 *
	 * It copy the first grid into the given grid (No ghost)
	 *
	 * \warning the Decomposition must be ensured to be the same, otherwise crashes can happen, if you want to copy the grid independently from the decomposition please use the operator equal
	 *
	 * \param g Grid to copy
	 *
	 * \return itself
	 *
	 */
	grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> & copy(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> & g)
	{
		auto it = this->getDomainIterator();

		while (it.isNext())
		{
			auto key = it.get();

			this->loc_grid.get(key.getSub()).get_o(key.getKey()) = g.loc_grid.get(key.getSub()).get_o(key.getKey());

			++it;
		}

		return *this;
	}

	/*! \brief Get the spacing on each dimension
	 *
	 * \return the spacing of the grid on each dimension as a point
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
	 * \param k grid_dist_key_dx point (in general returned by the iterators)
	 *
	 * \return the global position in the grid
	 *
	 */
	inline grid_key_dx<dim> getGKey(const grid_dist_key_dx<dim> & k)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
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
	 * \param output directory where to put the files + prefix
	 *
	 * \return true if the write operation succeed
	 *
	 */
	bool write(std::string output, size_t opt = VTK_WRITER | FORMAT_ASCII)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		file_type ft = file_type::ASCII;

		if (opt & FORMAT_BINARY)
			ft = file_type::BINARY;

		// Create a writer and write
		VTKWriter<boost::mpl::pair<device_grid,float>,VECTOR_GRIDS> vtk_g;
		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{
			Point<dim,St> offset = getOffset(i);
			vtk_g.add(loc_grid.get(i),offset,cd_sm.getCellBox().getP2(),gdb_ext.get(i).Dbox);
		}
		vtk_g.write(output + "_" + std::to_string(v_cl.getProcessUnitID()) + ".vtk", "grids", ft);

		return true;
	}

	/*! \brief Write the distributed grid information
	 *
	 * * grid_X.vtk Output each local grids for each local processor X
	 * * internal_ghost_X.vtk Internal ghost boxes in grid units for the local processor X
	 *
	 * \param output directory where to put the files + prefix
	 * \param i frame number
	 * \param opt options
	 *
	 * \return true id the write succeed
	 *
	 */
	bool write(std::string output, size_t i, size_t opt = VTK_WRITER | FORMAT_ASCII)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		file_type ft = file_type::ASCII;

		if (opt & FORMAT_BINARY)
			ft = file_type::BINARY;

		// Create a writer and write
		VTKWriter<boost::mpl::pair<device_grid,float>,VECTOR_GRIDS> vtk_g;
		for (size_t i = 0 ; i < loc_grid.size() ; i++)
		{
			Point<dim,St> offset = getOffset(i);
			vtk_g.add(loc_grid.get(i),offset,cd_sm.getCellBox().getP2(),gdb_ext.get(i).Dbox);
		}
		vtk_g.write(output + "_" + std::to_string(v_cl.getProcessUnitID()) + "_" + std::to_string(i) + ".vtk","grids",ft);

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


	/*! \brief It return the id of structure in the allocation list
	 *
	 * \see print_alloc and SE_CLASS2
	 *
	 * \return the id
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

	/*! \brief It print the internal ghost boxes and external ghost boxes in global unit
	 *
	 *
	 */
	void debugPrint()
	{
		std::cout << "-------- External Ghost boxes ---------- " << std::endl;

		for (size_t i = 0 ; i < eg_box.size() ; i++)
		{
			std::cout << "Processor: " << eg_box.get(i).prc << " Boxes:" << std::endl;

			for (size_t j = 0; j < eg_box.get(i).bid.size() ; j++)
			{
				std::cout << " Box: " << eg_box.get(i).bid.get(j).g_e_box.toString() << "   Id: " << eg_box.get(i).bid.get(j).g_id << std::endl;
			}
		}

		std::cout << "-------- Internal Ghost boxes ---------- " << std::endl;

		for (size_t i = 0 ; i < ig_box.size() ; i++)
		{
			std::cout << "Processor: " << ig_box.get(i).prc << " Boxes:" << std::endl;

			for (size_t j = 0 ; j < ig_box.get(i).bid.size() ; j++)
			{
				std::cout << " Box: " << ig_box.get(i).bid.get(j).box.toString() << "   Id: " << ig_box.get(i).bid.get(j).g_id << std::endl;
			}
		}
	}

	//! Define friend classes
	//\cond
	friend grid_dist_id<dim,St,T,typename Decomposition::extended_type,Memory,device_grid>;
	//\endcond
};



#endif
