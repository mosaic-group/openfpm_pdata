/*
 * grid_dist_id_iterator_sub.hpp
 *
 *  Created on: Feb 4, 2015
 *      Author: Pietro Incardona
 */

#ifndef GRID_DIST_ID_ITERATOR_HPP_
#define GRID_DIST_ID_ITERATOR_HPP_

#define FREE 1
#define FIXED 2
#define ITERATION_ISOLATION 4

#include "Grid/grid_dist_key.hpp"
#include "VCluster/VCluster.hpp"
#include "util/GBoxes.hpp"

#ifdef __NVCC__
#include "SparseGridGpu/encap_num.hpp"
#endif

template<unsigned int dim>
struct launch_insert_sparse_lambda_call
{
	template<typename ec_type, typename lambda_t,typename coord_type>
	__device__ inline static void call(ec_type & ec,lambda_t f, coord_type coord)
	{
		printf("Not implemented in this direction \n");
	}

	template<typename ite_type>
	__device__ inline static bool set_keys(grid_key_dx<3,int> & key, grid_key_dx<3,int> & keyg, ite_type & itg)
	{
		return false;
	}
};

template<>
struct launch_insert_sparse_lambda_call<3>
{
	template<typename grid_type, typename lambda_t1, typename lambda_t2,typename itd_type, typename coord_type>
	__device__ inline static void call(grid_type & grid,
									   lambda_t1 f1, lambda_t2 f2,
									   unsigned int blockId,
									   itd_type itd,
									   coord_type & key,
									   coord_type & keyg,unsigned int offset, bool & is_block_empty)
	{
#ifdef __NVCC__

	    bool is_active = f1(keyg.get(0),keyg.get(1),keyg.get(2));
	    is_active &= key.get(0) >= itd.start_base.get(0) && key.get(1) >= itd.start_base.get(1) && key.get(2) >= itd.start_base.get(2);

	    if (is_active == true)
	    {is_block_empty = false;}

	    __syncthreads();

	    if (is_block_empty == false)
	    {
	    	auto ec = grid.insertBlock(blockId);
	    	enc_num<decltype(grid.insertBlock(blockId))> ecn(ec,offset);

	        if ( is_active == true)
	        {
	        	f2(ecn,keyg.get(0),keyg.get(1),keyg.get(2));
	        	ec.template get<grid_type::pMask>()[offset] = 1;
	        }
	    }

#endif
	}

	template<typename ite_type>
	__device__ inline static bool set_keys(grid_key_dx<3,int> & key, grid_key_dx<3,int> & keyg, ite_type & itg)
	{
#ifdef __NVCC__

		key.set_d(0,threadIdx.x + blockIdx.x * blockDim.x + itg.start.get(0));
		key.set_d(1,threadIdx.y + blockIdx.y * blockDim.y + itg.start.get(1));
		key.set_d(2,threadIdx.z + blockIdx.z * blockDim.z + itg.start.get(2));

		keyg.set_d(0,key.get(0) + itg.origin.get(0));
		keyg.set_d(1,key.get(1) + itg.origin.get(1));
		keyg.set_d(2,key.get(2) + itg.origin.get(2));

		if (key.get(0) > itg.stop.get(0) || key.get(1) > itg.stop.get(1) || key.get(2) > itg.stop.get(2))
		{return true;}
#endif
		return false;
	}
};

template<>
struct launch_insert_sparse_lambda_call<2>
{
	template<typename grid_type, typename lambda_t1, typename lambda_t2,typename itd_type, typename coord_type>
	__device__ inline static void call(grid_type & grid,
									   lambda_t1 f1, lambda_t2 f2,
									   unsigned int blockId,
									   itd_type itd,
									   coord_type & key,
									   coord_type & keyg,unsigned int offset, bool & is_block_empty)
	{
#ifdef __NVCC__

	    bool is_active = f1(keyg.get(0),keyg.get(1));
	    is_active &= key.get(0) >= itd.start_base.get(0) && key.get(1) >= itd.start_base.get(1);

	    if (is_active == true)
	    {is_block_empty = false;}

	    __syncthreads();

	    if (is_block_empty == false)
	    {
	    	auto ec = grid.insertBlock(blockId);
	    	enc_num<decltype(grid.insertBlock(blockId))> ecn(ec,offset);

	        if ( is_active == true)
	        {
	        	f2(ecn,keyg.get(0),keyg.get(1));
	        	ec.template get<grid_type::pMask>()[offset] = 1;
	        }
	    }

#endif
	}

	template<typename ite_type>
	__device__ inline static bool set_keys(grid_key_dx<2,int> & key, grid_key_dx<2,int> & keyg, ite_type & itg)
	{
#ifdef __NVCC__
		key.set_d(0,threadIdx.x + blockIdx.x * blockDim.x + itg.start.get(0));
		key.set_d(1,threadIdx.y + blockIdx.y * blockDim.y + itg.start.get(1));

		keyg.set_d(0,key.get(0) + itg.origin.get(0));
		keyg.set_d(1,key.get(1) + itg.origin.get(1));

		if (key.get(0) > itg.stop.get(0) || key.get(1) > itg.stop.get(1))
		{return true;}
#endif
		return false;
	}
};

struct launch_insert_sparse
{
	template<typename grid_type, typename ite_type, typename lambda_f1, typename lambda_f2>
	__device__ void operator()(grid_type & grid, ite_type itg, bool & is_block_empty, lambda_f1 f1, lambda_f2 f2)
	{
#ifdef __NVCC__

		grid_key_dx<grid_type::dims,int> key;
		grid_key_dx<grid_type::dims,int> keyg;

		if (launch_insert_sparse_lambda_call<grid_type::dims>::set_keys(key,keyg,itg) == true)	{return;}

	    if (threadIdx.x == 0 && threadIdx.y == 0 && threadIdx.z == 0)
	    {is_block_empty = true;}

	    grid.init();

	    int offset = 0;
	    grid_key_dx<grid_type::dims,int> blk;
	    bool out = grid.template getInsertBlockOffset<ite_type>(itg,key,blk,offset);

	    auto blockId = grid.getBlockLinId(blk);

	    launch_insert_sparse_lambda_call<grid_type::dims>::call(grid,f1,f2,blockId,itg,key,keyg,offset,is_block_empty);

	    __syncthreads();

	    grid.flush_block_insert();
#endif
	}
};

template<bool is_free>
struct selvg
{
	template<typename a_it_type, typename gdb_ext_type, typename gList_type>
	static inline void call(a_it_type & a_it, gdb_ext_type & gdb_ext, gList_type & gList, size_t & g_c)
	{
		if (gdb_ext.get(g_c).Dbox.isValid() == false)
		{g_c++;}
		else
		{
			a_it.reinitialize(gList.get(g_c).getIterator(gdb_ext.get(g_c).Dbox.getKP1(),gdb_ext.get(g_c).Dbox.getKP2()));
			if (a_it.isNext() == false)	{g_c++;}
		}
	}
};

template<>
struct selvg<false>
{
	template<typename a_it_type, typename gdb_ext_type, typename gList_type>
	static inline void call(a_it_type & a_it, gdb_ext_type & gdb_ext, gList_type & gList, size_t & g_c)
	{
		// Full iterator (no subset)
		a_it.reinitialize(gList.get(g_c).getIterator());
		if (a_it.isNext() == false)	{g_c++;}
	}
};

/*! \brief Distributed grid iterator
 *
 * Iterator across the local elements of the distributed grid
 *
 * \tparam dim dimensionality of the grid
 * \tparam device_grid type of basic grid
 * \tparam stencil it inject the code to calculate stencil offset
 * \tparam sub_iterator it indicate the sub-iterator type of the device_grid
 *
 */
template<unsigned int dim, typename device_grid, typename device_sub_it, int impl, typename stencil  = no_stencil >
class grid_dist_iterator
{
	//! grid list counter
	size_t g_c;

	//! List of the grids we are going to iterate
	const openfpm::vector<device_grid> & gList;

	//! Extension of each grid: domain and ghost + domain
	const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext;

	//! Actual iterator
	device_sub_it a_it;

	//! stop point (is the grid size)
	grid_key_dx<dim> stop;

	/*! \brief from g_c increment g_c until you find a valid grid
	 *
	 */
	void selectValidGrid()
	{
		do
		{
			if (impl == FREE)
			{
				// When the grid has size 0 potentially all the other informations are garbage
				while (g_c < gList.size() && (gList.get(g_c).size() == 0 || gdb_ext.get(g_c).Dbox.isValid() == false ) ) g_c++;
			}
			else
			{
				// When the grid has size 0 potentially all the other informations are garbage
				while (g_c < gList.size() && (gList.get(g_c).size() == 0 || gdb_ext.get(g_c).GDbox.isValid() == false) ) g_c++;
			}

			// get the next grid iterator
			if (g_c < gList.size())
			{
				selvg<impl == FREE>::call(a_it,gdb_ext,gList,g_c);
			}
		} while (g_c < gList.size() && a_it.isNext() == false);

	}

	public:

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param gk std::vector of the local grid
	 * \param gdb_ext set of local subdomains
	 * \param stop end point
	 *
	 */
	grid_dist_iterator(const openfpm::vector<device_grid> & gk,
					   const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
					   const grid_key_dx<dim> & stop)
	:g_c(0),gList(gk),gdb_ext(gdb_ext),stop(stop)
	{
		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}


	/*! \brief Constructor of the distributed grid iterator with
	 *         stencil support
	 *
	 * \param gk std::vector of the local grid
	 * \param gdb_ext set of local subdomains
	 * \param stop end point
	 * \param stencil_pnt stencil points
	 *
	 */
	grid_dist_iterator(openfpm::vector<device_grid> & gk,
			           const openfpm::vector<GBoxes<device_grid::dims>> & gdb_ext,
					   const grid_key_dx<dim> & stop,
					   const grid_key_dx<dim> (& stencil_pnt)[stencil::nsp])
	:g_c(0),gList(gk),gdb_ext(gdb_ext),a_it(stencil_pnt),stop(stop)
	{
		// Initialize the current iterator
		// with the first grid
		selectValidGrid();
	}

	//! Copy constructor
	grid_dist_iterator(const grid_dist_iterator<dim,device_grid,device_sub_it,impl,stencil> & g)
	:g_c(g.g_c),gList(g.gList),gdb_ext(g.gdb_ext),a_it(g.a_it),stop(g.stop)
	{}

	//! Copy constructor
	grid_dist_iterator(grid_dist_iterator<dim,device_grid,device_sub_it,impl,stencil> && g)
	:g_c(g.g_c),gList(g.gList),gdb_ext(g.gdb_ext),a_it(g.a_it),stop(g.stop)
	{}

	//! Destructor
	~grid_dist_iterator()
	{
	}

	/*! \brief Get the next element
	 *
	 * \return the next grid_key
	 *
	 */
	inline grid_dist_iterator<dim,device_grid,device_sub_it,impl,stencil> & operator++()
	{
		++a_it;

		// check if a_it is at the end

		if (a_it.isNext() == true)
			return *this;
		else
		{
			// switch to the new grid
			g_c++;

			selectValidGrid();
		}

		return *this;
	}

	/*! \brief Check if there is the next element
	 *
	 * \return true if there is the next, false otherwise
	 *
	 */
	inline bool isNext() const
	{
		// If there are no other grid stop

		if (g_c >= gList.size())
		{return false;}

		return true;
	}

	/*! \brief Get the actual key
	 *
	 * \return the actual key
	 *
	 */
	inline grid_dist_key_dx<dim, typename device_grid::base_key> get() const
	{
		return grid_dist_key_dx<dim,typename device_grid::base_key>(g_c,a_it.get());
	}

	/*! \brief it return the stop point of the iterator
	 *
	 * The stop point of the iterator is just the grid size
	 *
	 * \return the stop point
	 *
	 */
	inline grid_key_dx<dim> getStop() const
	{
		return stop;
	}

	/*! \brief it return the start point of the iterator
	 *
	 * The start point of the iterator is the point with all coordinates zeros
	 *
	 * \return the start point
	 *
	 */
	inline grid_key_dx<dim> getStart() const
	{
		grid_key_dx<dim> start;

		start.zero();

		return start;
	}

	/*! \brief Get the boxes
	 *
	 *  Get the boxes that define the local grids
	 *
	 * \return Vector of local boxes
	 *
	 */
	inline const openfpm::vector<GBoxes<device_grid::dims>> & getGBoxes()
	{
		return gdb_ext;
	}

	/*! \brief Convert a g_dist_key_dx into a global key
	 *
	 * \see grid_dist_key_dx
	 * \see grid_dist_iterator
	 *
	 * \param k key position in local coordinates
	 *
	 * \return the global position in the grid
	 *
	 */
	inline grid_key_dx<dim> getGKey(const grid_dist_key_dx<dim,typename device_grid::base_key> & k)
	{
		// Get the sub-domain id
		size_t sub_id = k.getSub();

		auto k_glob = k.getKey();

		// shift
		auto k_glob2 = k_glob + gdb_ext.get(sub_id).origin;

		return k_glob2;
	}

	/*! \brief Return the stencil point offset
	 *
	 * \tparam id
	 *
	 * \return linearized distributed key
	 *
	 */
	template<unsigned int id> inline grid_dist_lin_dx getStencil()
	{
		return grid_dist_lin_dx(g_c,a_it.template getStencil<id>());
	}
};



#endif /* GRID_DIST_ID_ITERATOR_SUB_HPP_ */
