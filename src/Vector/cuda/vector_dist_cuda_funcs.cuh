/*
 * vector_dist_cuda_funcs.cuh
 *
 *  Created on: Aug 14, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_CUDA_FUNCS_CUH_
#define VECTOR_DIST_CUDA_FUNCS_CUH_

#include "Vector/util/vector_dist_funcs.hpp"
#include "util/cuda/moderngpu/kernel_reduce.hxx"

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

template<unsigned int dim, typename St, typename cartdec_gpu, typename particles_type, typename vector_out>
__global__ void process_id_proc_each_part(cartdec_gpu cdg, particles_type parts, vector_out output , int rank)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    cdg.applyPointBC(parts.get(p));
	Point<dim,St> xp = parts.template get<0>(p);

	int pr = cdg.processorID(xp);

	output.template get<1>(p) = (pr == rank)?-1:pr;
	output.template get<0>(p) = p;
}

template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= vd.size() - 1) return;

    if (vd.template get<prp_off>(p) != vd.template get<prp_off>(p+1))
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

template<typename vector_g_opart_type, typename vector_prp_type_out, typename vector_prp_type_in, unsigned int ... prp>
__global__  void process_ghost_particles_prp(vector_g_opart_type g_opart, vector_prp_type_out m_prp,
		     	 	 	 	 	 	   	     vector_prp_type_in  v_prp, unsigned int offset)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= m_prp.size()) return;

    process_ghost_device_particle_prp<vector_g_opart_type,vector_prp_type_out,vector_prp_type_in,prp...>(i,offset,g_opart,m_prp,v_prp);
}

template<unsigned int dim, typename vector_g_opart_type, typename vector_pos_type_out, typename vector_pos_type_in, typename vector_shift_type_in>
__global__  void process_ghost_particles_pos(vector_g_opart_type g_opart, vector_pos_type_out m_pos,
		     	 	 	 	 	 	   	     vector_pos_type_in  v_pos, vector_shift_type_in shifts, unsigned int offset)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= m_pos.size()) return;

    unsigned long int psid = g_opart.template get<1>(i+offset);

    unsigned int id = psid & 0xFFFFFFFF;
    unsigned int shift_id = psid >> 32;

	#pragma unroll
    for (int j = 0; j < dim ; j++)
    {
    	m_pos.template get<0>(i)[j] = v_pos.template get<0>(id)[j] - shifts.template get<0>(shift_id)[j];
    }
}

template<bool with_pos, unsigned int dim, typename vector_g_opart_type, typename vector_pos_type,
         typename vector_prp_type, typename vector_shift_type_in>
__global__ void process_ghost_particles_local(vector_g_opart_type g_opart, vector_pos_type v_pos, vector_prp_type v_prp,
											  vector_shift_type_in shifts, unsigned int base)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= g_opart.size()) return;

    unsigned int pid = g_opart.template get<0>(i);
    unsigned int shift_id = g_opart.template get<1>(i);

    if (with_pos == true)
    {
		#pragma unroll
		for (int j = 0; j < dim ; j++)
		{
			v_pos.template get<0>(base+i)[j] = v_pos.template get<0>(pid)[j] - shifts.template get<0>(shift_id)[j];
		}
    }

    v_prp.set(base+i,v_prp.get(pid));
}

template<unsigned int dim, typename St, typename vector_of_box, typename vector_type,  typename output_type>
__global__ void num_shift_ghost_each_part(vector_of_box box_f, vector_type vd,  output_type out)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= vd.size()) return;

    Point<dim,St> xp = vd.template get<0>(p);

    unsigned int n = 0;

    for (unsigned int i = 0 ; i < box_f.size() ; i++)
    {
    	if (Box<dim,St>(box_f.get(i)).isInsideNP(xp) == true)
    	{n++;}
    }

    out.template get<0>(p) = n;
}

template<unsigned int dim, typename St,
         typename vector_of_box,
         typename vector_of_shifts,
         typename vector_type_pos,
         typename vector_type_prp,
         typename start_type,
         typename shifts_type,
         typename output_type>
__global__ void shift_ghost_each_part(vector_of_box box_f, vector_of_shifts box_f_sv,
		                              vector_type_pos v_pos, vector_type_prp v_prp,
		                              start_type start, shifts_type shifts,
		                              output_type output, unsigned int offset)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= v_pos.size()) return;

    Point<dim,St> xp = v_pos.template get<0>(p);

    unsigned int base_o = start.template get<0>(p);
    unsigned int base = base_o + offset;


    unsigned int n = 0;

    for (unsigned int i = 0 ; i < box_f.size() ; i++)
    {
    	if (Box<dim,St>(box_f.get(i)).isInsideNP(xp) == true)
    	{
    		unsigned int shift_id = box_f_sv.template get<0>(i);

#pragma unroll
    		for (unsigned int j = 0 ; j < dim ; j++)
    		{
    			v_pos.template get<0>(base+n)[j] = xp.get(j) - shifts.template get<0>(shift_id)[j];
    			output.template get<0>(base_o+n) = p;
    			output.template get<1>(base_o+n) = shift_id;
    		}

    		v_prp.set(base+n,v_prp.get(p));

    		n++;
    	}
    }
}

template<unsigned int prp, typename vector_type>
auto reduce(vector_type & vd) -> typename std::remove_reference<decltype(vd.template getProp<prp>(0))>::type
{
	typedef typename std::remove_reference<decltype(vd.template getProp<prp>(0))>::type reduce_type;

	CudaMemory mem;
	mem.allocate(sizeof(reduce_type));

	mgpu::reduce((reduce_type *)vd.getPropVector(). template getDeviceBuffer<prp>(),
			            vd.size_local(), (reduce_type *)mem.getDevicePointer() ,
			            mgpu::plus_t<reduce_type>(), vd.getVC().getmgpuContext());

	mem.deviceToHost();

	return *(reduce_type *)(mem.getPointer());
}

#endif /* VECTOR_DIST_CUDA_FUNCS_CUH_ */
