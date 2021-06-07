/*
 * grid_dist_id_kernels.cuh
 *
 *  Created on: Jun 25, 2019
 *      Author: i-bird
 */

#ifndef GRID_DIST_ID_KERNELS_CUH_
#define GRID_DIST_ID_KERNELS_CUH_

#include "config.h"

#ifdef CUDA_GPU

template<unsigned int dim>
struct ite_gpu_dist
{
	dim3 thr;
	dim3 wthr;

	grid_key_dx<dim,int> start;
	grid_key_dx<dim,int> stop;

	grid_key_dx<dim,int> start_base;

	grid_key_dx<dim,int> origin;

	ite_gpu_dist()
	{}

	ite_gpu_dist(ite_gpu<dim> & ite)
	{
		thr = ite.thr;
		wthr = ite.wthr;

		start = ite.start;
		stop = ite.stop;
	}

	size_t nblocks()
	{
		return wthr.x * wthr.y * wthr.z;
	}
};

#define GRID_ID_3_GLOBAL(ite_gpu) grid_key_dx<3,int> key;\
								  grid_key_dx<3,int> keyg;\
							  key.set_d(0,threadIdx.x + blockIdx.x * blockDim.x + ite_gpu.start.get(0));\
    						  key.set_d(1,threadIdx.y + blockIdx.y * blockDim.y + ite_gpu.start.get(1));\
							  key.set_d(2,threadIdx.z + blockIdx.z * blockDim.z + ite_gpu.start.get(2));\
							  \
							  bool inactive = false;\
							  \
							  keyg.set_d(0,key.get(0) + ite_gpu.origin.get(0));\
    						  keyg.set_d(1,key.get(1) + ite_gpu.origin.get(1));\
							  keyg.set_d(2,key.get(2) + ite_gpu.origin.get(2));\
										 \
										 if (key.get(0) > ite_gpu.stop.get(0) || key.get(1) > ite_gpu.stop.get(1) || key.get(2) > ite_gpu.stop.get(2))\
    									 {inactive = true;}


#define GRID_ID_2_GLOBAL(ite_gpu) grid_key_dx<2,int> key;\
								  grid_key_dx<2,int> keyg;\
							  key.set_d(0,threadIdx.x + blockIdx.x * blockDim.x + ite_gpu.start.get(0));\
    						  key.set_d(1,threadIdx.y + blockIdx.y * blockDim.y + ite_gpu.start.get(1));\
							  \
							  bool inactive = false;\
							  \
							  keyg.set_d(0,key.get(0) + ite_gpu.origin.get(0));\
    						  keyg.set_d(1,key.get(1) + ite_gpu.origin.get(1));\
										 \
										 if (key.get(0) > ite_gpu.stop.get(0) || key.get(1) > ite_gpu.stop.get(1))\
    									 {inactive = true;}

#endif

template<typename grid_type, typename ite_gpu_type,typename func_t,typename ... args_t>
__global__ void grid_apply_functor(grid_type g, ite_gpu_type ite, func_t f, args_t ... args)
{
	f(g,ite,args...);
}

template<typename grid_type, typename ite_gpu_type,typename func_t,typename ... args_t>
__global__ void grid_apply_functor_shared_bool(grid_type g, ite_gpu_type ite, func_t f, args_t ... args)
{
	__shared__ bool is_empty_block;

	f(g,ite,is_empty_block,args...);
}

#endif /* GRID_DIST_ID_KERNELS_CUH_ */
