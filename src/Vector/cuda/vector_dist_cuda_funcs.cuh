/*
 * vector_dist_cuda_funcs.cuh
 *
 *  Created on: Aug 14, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_CUDA_FUNCS_CUH_
#define VECTOR_DIST_CUDA_FUNCS_CUH_

#include "Vector/util/vector_dist_funcs.hpp"

template<typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= vd.size() - 1) return;

    if (vd.template get<0>(p) != vd.template get<0>(p+1))
	{
    	int i = atomicAdd(cnt, 1);
    	offs.template get<0>(i) = p+1;
    	offs.template get<1>(i) = vd.template get<0>(p);
	}
}

template<typename vector_m_opart_type, typename vector_pos_type_out, typename vector_prp_type_out,
		 typename vector_pos_type_in,  typename vector_prp_type_in>
__global__  void process_map_particles(vector_m_opart_type m_opart, vector_pos_type_out m_pos, vector_prp_type_out m_prp,
		     	 	 	 	 	 	   vector_pos_type_in  v_pos, vector_prp_type_in  v_prp, unsigned int offset)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= m_pos.size()) return;

    process_map_device_particle<proc_without_prp_device>(i,offset,m_opart,m_pos,m_prp,v_pos,v_prp);
}


#endif /* VECTOR_DIST_CUDA_FUNCS_CUH_ */
