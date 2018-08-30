/*
 * vector_dist_cuda_funcs.cuh
 *
 *  Created on: Aug 14, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_CUDA_FUNCS_CUH_
#define VECTOR_DIST_CUDA_FUNCS_CUH_

#include "Vector/util/vector_dist_funcs.hpp"

template<unsigned int dim, typename St, typename decomposition_type, typename vector_type, typename start_type, typename output_type>
__global__ void proc_label_id_ghost(decomposition_type dec,vector_type vd, start_type starts, output_type out)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= vd.size()) return;

    Point<dim,St> xp = vd.template get<0>(p);

    unsigned int base = starts.template get<0>(p);

    dec.ghost_processor_ID(xp,out,base,p);
}

template<unsigned int dim, typename St, typename decomposition_type, typename vector_type,  typename output_type>
__global__ void num_proc_ghost_each_part(decomposition_type dec, vector_type vd,  output_type out)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= vd.size()) return;

    Point<dim,St> xp = vd.template get<0>(p);

    out.template get<0>(p) = dec.ghost_processorID_N(xp);
}

template<typename cartdec_gpu, typename particles_type, typename vector_out>
__global__ void process_id_proc_each_part(cartdec_gpu cdg, particles_type parts, vector_out output , int rank)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

	Point<3,float> xp = parts.template get<0>(p);

	int pr = cdg.processorIDBC(xp);

	output.template get<1>(p) = (pr == rank)?-1:pr;
	output.template get<0>(p) = p;
}

template<typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= vd.size() - 1) return;

    if (vd.template get<1>(p) != vd.template get<1>(p+1))
	{
    	int i = atomicAdd(cnt, 1);
    	offs.template get<0>(i) = p+1;
    	offs.template get<1>(i) = vd.template get<1>(p);
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
