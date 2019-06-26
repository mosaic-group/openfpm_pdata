/*
 * grid_dist_id_kernels.cuh
 *
 *  Created on: Jun 25, 2019
 *      Author: i-bird
 */

#ifndef GRID_DIST_ID_KERNELS_CUH_
#define GRID_DIST_ID_KERNELS_CUH_


template<typename grid_type, typename func_t,typename ... args_t>
__global__ void grid_apply_functor(grid_type g, ite_gpu<grid_type::d> ite, func_t f, args_t ... args)
{
	f(g,ite,args...);
}



#endif /* GRID_DIST_ID_KERNELS_CUH_ */
