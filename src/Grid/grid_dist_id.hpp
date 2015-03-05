#ifndef COM_UNIT_HPP
#define COM_UNIT_HPP

#include <vector>
#include "Grid/map_grid.hpp"
#include "VCluster.hpp"
#include "Space/SpaceBox.hpp"
#include "mathutil.hpp"
#include "grid_dist_id_iterator.hpp"
#include "grid_dist_key.hpp"

#define SUB_UNIT_FACTOR 64


/*! \brief This is a distributed grid
 *
 * Implementation of a distributed grid with id decomposition. A distributed grid is a grid distributed
 * across processors. The decomposition is performed on the id of the elements
 *
 * [Examples]
 *
 * on 1D where the id is from 1 to N
 * processor k take M contiguous elements
 *
 * on 3D where (for example)
 * processor k take M id-connected elements
 *
 * \param dim Dimensionality of the grid
 * \param T type of grid
 * \param Decomposition Class that decompose the grid for example CartDecomposition
 * \param Mem Is the allocator
 * \param device type of base structure is going to store the data
 *
 */

template<unsigned int dim, typename T, typename Decomposition,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_id
{
	// Ghost expansion
	Box<dim,size_t> ghost;

	//! Local grids
	Vcluster_object_array<device_grid> loc_grid;

	//! Space Decomposition
	Decomposition dec;

	//! Size of the grid on each dimension
	size_t g_sz[dim];

	//! Communicator class

	Vcluster & v_cl;

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
			v_size[d] = sp.getHigh(d) - sp.getLow(d) + 1;
		}
	}

public:

	//! constructor
	grid_dist_id(Vcluster v_cl, Decomposition & dec, size_t (& g_sz)[dim], Box<dim,size_t> & ghost)
	:ghost(ghost),loc_grid(NULL),v_cl(v_cl),dec(dec)
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
		{div[i] = round_big_2(pow(n_sub,1.0/dim));}

		// Create the sub-domains
		dec.setParameters(div);

		// Create local grid
		Create();
	}

	//! constructor
	grid_dist_id(size_t (& g_sz)[dim])
	:dec(Decomposition(*global_v_cluster)),v_cl(*global_v_cluster)
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
		{div[i] = round_big_2(pow(n_sub,1.0/dim));}

		// Box
		Box<dim,size_t> b(g_sz);

		// Create the sub-domains
		dec.setParameters(div,b);

		// Create local grid
		Create();
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
			// Get the local hyper-cube

			SpaceBox<dim,size_t> sp = dec.getLocalHyperCube(i);

			// Calculate the local grid size

			getGridSize(sp,l_res);

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
		grid_dist_iterator<dim,device_grid> it(loc_grid,0);

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
};

#endif
