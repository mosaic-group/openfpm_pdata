/*
 * grid_dist_id_iterator_dec_gpu.cuh
 *
 *  Created on: Sep 1, 2019
 *      Author: i-bird
 */

#ifndef GRID_DIST_ID_ITERATOR_DEC_GPU_CUH_
#define GRID_DIST_ID_ITERATOR_DEC_GPU_CUH_


#include "Grid/Iterators/grid_dist_id_iterator.hpp"
#include "Grid/grid_dist_util.hpp"
#include "Grid/Iterators/grid_dist_id_iterator_util.hpp"
#include "Grid/cuda/grid_dist_id_kernels.cuh"

template<unsigned int impl>
struct launch_call_impl
{
	template<typename loc_grid_type, typename ite_type, typename itd_type, typename functor_type, typename ... argsT>
	inline static void call(loc_grid_type & loc_grid, ite_type & ite , itd_type & itd, functor_type functor, argsT ... args)
	{
		CUDA_LAUNCH(grid_apply_functor,ite,loc_grid.toKernel(), itd, functor, args... );
	}
};

template<>
struct launch_call_impl<1>
{
	template<typename loc_grid_type, typename ite_type, typename itd_type, typename functor_type,typename ... argsT>
	inline static void call(loc_grid_type & loc_grid, ite_type & ite, itd_type & itd, functor_type f, argsT ... args)
	{

		auto g = loc_grid.toKernel();

		auto lamb = [g,itd,f,args ...] __device__ () mutable
		{
			__shared__ bool is_empty_block;

			f(g,itd,is_empty_block,args...);
		};

		CUDA_LAUNCH_LAMBDA_TLS(ite,lamb);
	}
};

/*! \brief Given the decomposition it create an iterator
 *
 * Iterator across the local elements of the distributed grid
 *
 * \tparam dec Decomposition type
 *
 */
template<typename Decomposition, typename deviceGrids, bool ghost_or_domain = false>
class grid_dist_id_iterator_gpu
{
	//! grid list counter
	size_t g_c;

	//! Extension of each grid: domain and ghost + domain
	openfpm::vector<GBoxes<Decomposition::dims>> gdb_ext;

	//! start key
	grid_key_dx<Decomposition::dims> start;

	//! stop key
	grid_key_dx<Decomposition::dims> stop;

	//! Local device grids
	deviceGrids & loc_grids;

	//! number of threads to launch the kernels
	size_t n_thr;

	//! Maximum number of insertions for each GPU block
	int nSlot = -1;

	//! Spacing
	typename Decomposition::stype spacing[Decomposition::dims];

	public:

	/*! \brief Copy operator=
	*
	* \param tmp iterator to copy
	*
	* \return itself
	*
	*/
//	grid_dist_id_iterator_gpu<Decomposition,deviceGrids> & operator=(const grid_dist_id_iterator_gpu<Decomposition,deviceGrids> & tmp)
//	{
//		g_c = tmp.g_c;
//		gdb_ext = tmp.gdb_ext;

//		start = tmp.start;
//		stop = tmp.stop;
//		loc_grids = tmp.loc_grids;

//		return *this;
//	}

	/*! \brief Copy constructor
	*
	* \param tmp iterator to copy
	*
	*/
//	grid_dist_id_iterator_gpu(const grid_dist_id_iterator_gpu<Decomposition,deviceGrids> & tmp)
//	:loc_grids(tmp.loc_grids)
//	{
//		this->operator=(tmp);
//	}

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param dec Decomposition
	 * \param sz size of the grid
	 *
	 */
	grid_dist_id_iterator_gpu(deviceGrids & loc_grids,Decomposition & dec, const size_t (& sz)[Decomposition::dims])
	:loc_grids(loc_grids),g_c(0)
	{
		// Initialize start and stop
		start.zero();
		for (size_t i = 0 ; i < Decomposition::dims ; i++)
			stop.set_d(i,sz[i]-1);

		// From the decomposition construct gdb_ext
		create_gdb_ext<Decomposition::dims,Decomposition>(gdb_ext,dec,sz,dec.getDomain(),spacing);

		g_c = 0;
	}

	/*! \brief Constructor of the distributed grid iterator
	 *
	 * \param dec Decomposition
	 * \param sz size of the grid
	 * \param start point
	 * \param stop point
	 *
	 */
	grid_dist_id_iterator_gpu(deviceGrids & loc_grids ,Decomposition & dec, const size_t (& sz)[Decomposition::dims], grid_key_dx<Decomposition::dims> start, grid_key_dx<Decomposition::dims> stop)
	:loc_grids(loc_grids),g_c(0),start(start),stop(stop)
	{
		// From the decomposition construct gdb_ext
		create_gdb_ext<Decomposition::dims,Decomposition>(gdb_ext,dec,sz,dec.getDomain(),spacing);

		g_c = 0;
	}

	// Destructor
	~grid_dist_id_iterator_gpu()
	{
	}

	/*! \brief The the number of maximum inserts each GPU block can do
	 *
	 * \param nSlot maximum number of insert
	 *
	 */
	void setGPUInsertBuffer(int nSlot)
	{
		this->nSlot = nSlot;
	}

	/*! \brief Set the number of threads for each block
	 *
	 * \param nthr number of threads for each block
	 *
	 */
	void setBlockThreads(size_t nthr)
	{
		this->n_thr = nthr;
	}

	/*! \brief Return true if we point to a valid grid
	 *
	 * \return true if valid grid
	 *
	 */
	inline bool isNextGrid()
	{
		return g_c < gdb_ext.size();
	}

	/*! \brief Return the index of the grid in which we are iterating
	 *
	 *
	 */
	inline size_t getGridId()
	{
		return g_c;
	}

	/*! \brief next grid
	 *
	 *
	 */
	inline void nextGrid()
	{
		g_c++;
	}


	/*! \brief Get the spacing of the grid
	 *
	 * \param i
	 *
	 */
	inline typename Decomposition::stype getSpacing(size_t i)
	{
		return spacing[i];
	}

	/*! \brief Launch a functor with a particular kernel
	 *
	 * \param functor function kernel
	 * \param argsType arguments
	 *
	 */
	template<unsigned int impl = 0, typename func_t, typename ... argsType >
	inline void launch(func_t functor,argsType ... args)
	{
		for (g_c = 0 ; g_c < gdb_ext.size() ; g_c++)
		{
			ite_gpu_dist<Decomposition::dims> itd;
			ite_gpu<Decomposition::dims> ite;

			// intersect

			Box<Decomposition::dims,int> range_box(start,stop);
			Box<Decomposition::dims,int> kbox;
			range_box -= gdb_ext.get(g_c).origin;
			bool intersect = range_box.Intersect(gdb_ext.get(g_c).Dbox,kbox);

			if (intersect == false)	{continue;}

			auto & lg = loc_grids.get(g_c);

			for (int i = 0 ; i < Decomposition::dims ; i++)
			{
				ite.start.set_d(i,(kbox.getKP1().get(i) / lg.getBlockEdgeSize())*lg.getBlockEdgeSize() );
				ite.stop.set_d(i,  kbox.getKP2().get(i));
			}

			// the thread extensions are

			for (int i = 0 ; i < Decomposition::dims ; i++)
			{
				itd.origin.set_d(i,gdb_ext.get(g_c).origin.get(i));
				itd.start_base.set_d(i,kbox.getKP1().get(i) % lg.getBlockEdgeSize() + ite.start.get(i));
			}

			ite.thr.x = lg.getBlockEdgeSize();
			ite.wthr.x = (ite.stop.get(0) - ite.start.get(0) + 1) / lg.getBlockEdgeSize() + ((ite.stop.get(0) - ite.start.get(0) + 1) % lg.getBlockEdgeSize() != 0);

			ite.thr.y = lg.getBlockEdgeSize();
			ite.wthr.y = (ite.stop.get(1) - ite.start.get(1) + 1) / lg.getBlockEdgeSize() + ((ite.stop.get(1) - ite.start.get(1) + 1) % lg.getBlockEdgeSize() != 0);

			if (Decomposition::dims > 2)
			{
				ite.thr.z = lg.getBlockEdgeSize();
				ite.wthr.z = (ite.stop.get(2) - ite.start.get(2) + 1) / lg.getBlockEdgeSize() + ((ite.stop.get(2) - ite.start.get(2) + 1) % lg.getBlockEdgeSize() != 0);
			}

			itd.wthr = ite.wthr;
			itd.thr = ite.thr;
			itd.start = ite.start;
			itd.stop = ite.stop;

			if (nSlot != -1)
			{
				loc_grids.get(g_c).setGPUInsertBuffer((unsigned int)ite.nblocks(),(unsigned int)nSlot);
			}

			if (ite.nblocks() != 0)
			{
				launch_call_impl<impl>::call(loc_grids.get(g_c),ite,itd,functor,args...);
			}
		}
	}


	/*! \brief Get the starting point of the sub-grid we are iterating
	 *
	 * \return the starting point
	 *
	 */
	inline grid_key_dx<Decomposition::dims> getStart()
	{
		return start;
	}

	/*! \brief Get the starting point of the sub-grid we are iterating
	 *
	 * \return the stop point
	 *
	 */
	inline grid_key_dx<Decomposition::dims> getStop()
	{
		return stop;
	}
};


#endif /* GRID_DIST_ID_ITERATOR_DEC_GPU_CUH_ */
