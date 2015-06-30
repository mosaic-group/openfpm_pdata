#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include "Grid/map_grid.hpp"
#include "VCluster.hpp"
#include "Space/SpaceBox.hpp"
#include "mathutil.hpp"
#include "grid_dist_id_iterator.hpp"
#include "grid_dist_key.hpp"
#include "NN/CellList/CellDecomposer.hpp"
#include "util/object_util.hpp"
#include "memory/ExtPreAlloc.hpp"
#include "VTKWriter.hpp"

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

	/*! \brief Create per-processor internal ghost box list in grid units
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

				Box_id bid_t(cvt);
				bid_t.id = 0/*g.LinId(bid_t.box.getKP1())*/;
				pib.bid.add(bid_t);
			}
		}

		init_i_g_box = true;
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

			// Center the local grid to zero
			sp_t -= sp_t.getP1();

			// save the domain box seen inside the domain + ghost box (see GDBoxes for a visual meaning)
			gdb_ext.last().Dbox = sp_t;
			gdb_ext.last().Dbox -= g_int_t.getP1();
			gdb_ext.last().Dbox.shrinkP2(1);

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

	/*! \brief it store a box and its unique id
	 *
	 */
	struct Box_id
	{
		//! Constructor
		inline Box_id(const ::Box<dim,size_t> & box)
		:box(box)
		{}

		//! Box
		::Box<dim,size_t> box;

		//! id
		size_t id;
	};

	/*! \brief Internal ghost box
	 *
	 */
	struct p_box_grid
	{
		// Internal ghost in grid units
		openfpm::vector<Box_id> bid;

		//! processor id
		size_t prc;
	};

	//! Flag that indicate if internal ghost box has been initialized
	bool init_i_g_box = false;

	//! Internal ghost boxes in grid units
	openfpm::vector<p_box_grid> ig_box;

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

		// send vector for each processor
		typedef  openfpm::vector<prp_object,openfpm::device_cpu<prp_object>,ExtPreAlloc<Memory>> send_vector;

		// Send buffer size in byte ( one buffer for all processors )
		size_t size_byte_prp = 0;

		create_ig_box();

		// Convert the ghost boxes into grid unit boxes
/*		const openfpm::vector<p_box_grid> p_box_g;

		// total number of sending vector
		size_t n_send_vector = 0;

		// Calculate the total size required for the sending buffer for each processor
		for ( size_t i = 0 ; i < p_box_g.size() ; i++ )
		{
			// for each ghost box
			for (size_t j = 0 ; j < p_box_g.get(i).box.size() ; j++)
			{
				size_t alloc_ele = openfpm::vector<prp_object>::calculateMem(p_box_g.get(i).box.get(j).volume(),0);
				pap_prp.push_back(alloc_ele);
				size_byte_prp += alloc_ele;

				n_send_vector++;
			}
		}

		// resize the property buffer memory
		g_prp_mem.resize(size_byte_prp);

		// Create an object of preallocated memory for properties
		ExtPreAlloc<Memory> * prAlloc_prp = new ExtPreAlloc<Memory>(pap_prp,g_prp_mem);

		// create a vector of send vector (ExtPreAlloc warrant that all the created vector are contiguous)
		openfpm::vector<send_vector> g_send_prp;

		// create a number of send buffers equal to the near processors
		g_send_prp.resize(n_send_vector);

		// Create the vectors giving them memory
		for ( size_t i = 0 ; i < p_box_g.size() ; i++ )
		{
			// for each ghost box
			for (size_t j = 0 ; j < p_box_g.get(i).box.size() ; j++)
			{
				// set the preallocated memory to ensure contiguity
				g_send_prp.get(i).setMemory(*prAlloc_prp);

				// resize the sending vector (No allocation is produced)
				g_send_prp.get(i).resize(ghost_prc_sz.get(i));
			}
		}

		// Fill the sending buffers and produce a sending request
		for ( size_t i = 0 ; i < p_box_g.size() ; i++ )
		{
			// for each ghost box
			for (size_t j = 0 ; j < p_box_g.get(i).box.size() ; j++)
			{
				// Create a sub grid iterator of the ghost
				grid_key_dx_iterator_sub<dim> g_it(,p_box_g.getKP1(),p_box_g.getKP2());

				while (g_it.isNext())
				{
					// copy all the object in the send buffer
					typedef encapc<1,prop,typename openfpm::vector<prop>::memory_t> encap_src;
					// destination object type
					typedef encapc<1,prp_object,typename openfpm::vector<prp_object>::memory_t> encap_dst;

					// Copy only the selected properties
					object_si_d<encap_src,encap_dst,ENCAP,prp...>(v_prp.get(INTERNAL).get(opart.get(i).get(j)),g_send_prp.get(i).get(j));

					++g_it;
				}
			}
			// Queue a send request
			v_cl.send(p_box_g.get(i).proc,g_send_prp.get(i));
		}

		// Calculate the receive pattern and allocate the receive buffers

		// For each processor get the ghost boxes
		const openfpm::vector<p_box> & p_box_ext = dec.getGhostExternalBox();

		// Convert the ghost boxes into grid unit boxes
		const openfpm::vector<p_box_grid> p_box_g_ext;

		size_byte_prp = 0;

		// Calculate the receive pattern
		for ( size_t i = 0 ; i < p_box_g_ext.size() ; i++ )
		{
			// for each ghost box
			for (size_t j = 0 ; j < p_box_g_ext.get(i).box.size() ; j++)
			{
				size_t alloc_ele = openfpm::vector<prp_object>::calculateMem(p_box_g_ext.get(i).box.get(j).volume(),0);
				pap_prp_recv.push_back(alloc_ele);
				size_byte_prp += alloc_ele;
			}

			v_cl.recv();
		}

		// Fill the internal ghost

		// wait to receive communication
		v_cl.execute();

		// move the received buffers into the ghost part

		 */
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
