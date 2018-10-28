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
#include "util/cuda/moderngpu/kernel_scan.hxx"
#include "Decomposition/common.hpp"

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

template<unsigned int dim, typename St, typename particles_type>
__global__ void apply_bc_each_part(Box<dim,St> domain, periodicity_int<dim> bc, particles_type parts)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    applyPointBC_no_dec(domain,bc,parts.get(p));
}

template<bool merge_pos, typename vector_pos_type, typename vector_prp_type, typename stns_type, unsigned int ... prp>
__global__ void merge_sort_part(vector_pos_type vd_pos, vector_prp_type vd_prp,
		                        vector_pos_type v_pos_ord, vector_prp_type vd_prp_ord,
		                        stns_type nss)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd_pos.size()) return;

	if (merge_pos == true)
	{
		vd_pos.template set<0>(p,v_pos_ord,nss.template get<0>(p));
	}

	vd_prp.template set<prp...>(p,vd_prp_ord,nss.template get<0>(p));
}

template<unsigned int dim, typename St, typename cartdec_gpu, typename particles_type, typename vector_out, typename prc_sz_type>
__global__ void process_id_proc_each_part(cartdec_gpu cdg, particles_type parts, vector_out output, prc_sz_type prc_sz , int rank)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= parts.size()) return;

    cdg.applyPointBC(parts.get(p));
	Point<dim,St> xp = parts.template get<0>(p);

	int pr = cdg.processorID(xp);

#ifndef TEST1
	output.template get<1>(p) = (pr == rank)?-1:pr;
	output.template get<0>(p) = p;
#else
	output.template get<1>(p) = pr;
	int nl = atomicAdd(&prc_sz.template get<0>(pr), 1);
	output.template get<2>(p) = nl;
#endif
}

template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)vd.size() - 1) return;

    if (vd.template get<prp_off>(p) != vd.template get<prp_off>(p+1))
	{
    	int i = atomicAdd(cnt, 1);
    	offs.template get<0>(i) = p+1;
    	offs.template get<1>(i) = vd.template get<1>(p);
	}
}

template<unsigned int prp_off, typename vector_type,typename vector_type_offs>
__global__  void find_buffer_offsets_no_prc(vector_type vd, int * cnt, vector_type_offs offs)
{
    int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= (int)vd.size() - 1) return;

    if (vd.template get<prp_off>(p) != vd.template get<prp_off>(p+1))
	{
    	int i = atomicAdd(cnt, 1);
    	offs.template get<0>(i) = p+1;
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

template<typename vector_lbl_type, typename starts_type>
__global__  void reorder_lbl(vector_lbl_type m_opart, starts_type starts)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= m_opart.size()) return;

    int pr = m_opart.template get<1>(i);

    m_opart.template get<0>(starts.template get<0>(pr) + m_opart.template get<2>(i)) = i;
}

template<typename red_type>
struct _add_: mgpu::plus_t<red_type>
{};

template<typename red_type>
struct _max_: mgpu::maximum_t<red_type>
{};

template<unsigned int prp, template <typename> class op, typename vector_type>
auto reduce(vector_type & vd) -> typename std::remove_reference<decltype(vd.template getProp<prp>(0))>::type
{
	typedef typename std::remove_reference<decltype(vd.template getProp<prp>(0))>::type reduce_type;

	CudaMemory mem;
	mem.allocate(sizeof(reduce_type));

	mgpu::reduce((reduce_type *)vd.getPropVector(). template getDeviceBuffer<prp>(),
			            vd.size_local(), (reduce_type *)mem.getDevicePointer() ,
			            op<reduce_type>(), vd.getVC().getmgpuContext());

	mem.deviceToHost();

	return *(reduce_type *)(mem.getPointer());
}

template<typename vector_type>
__global__ void create_index(vector_type vd)
{
	int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= vd.size()) return;

    vd.template get<0>(i) = i;
}

template<unsigned int dim, typename vector_pos_type, typename vector_prp_type, typename scan_type>
__global__ void copy_new_to_old(vector_pos_type vd_pos_dst, vector_prp_type vd_prp_dst, vector_pos_type vd_pos_src, vector_prp_type vd_prp_src, scan_type idx)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= vd_prp_dst.size()) return;

    for (unsigned int k = 0 ; k < dim ; k++)
    {vd_pos_dst.template get<0>(i)[k] = vd_pos_src.template get<0>(idx.template get<0>(i))[k];}

    vd_prp_dst.set(i,vd_prp_src,idx.template get<0>(i));
}

/*! \brief Remove the particles marked on the properties prp (particles marked has has property set to 1, the others to 0)
 *
 * \warning the function is destructive on prp, it mean that after destruction the prp of the particles can contain garbage
 *
 * \tparam prp property that indicate the particles to remove
 *
 * \param vd distributed vector
 *
 */
template<unsigned int prp, typename vector_type>
void remove_marked(vector_type & vd)
{
	// This function make sense only if prp is an int or unsigned int
	if (std::is_same< typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type, int >::value == false &&
		std::is_same< typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type, unsigned int >::value == false &&
		std::is_same< typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type, float >::value == false &&
		std::is_same< typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type, double >::value == false &&
		std::is_same< typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type, size_t >::value == false)
	{
		std::cout << __FILE__ << ":" << __LINE__ << " error, the function remove_marked work only if is an integer or unsigned int" << std::endl;
		return;
	}

	typedef typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type remove_type;

	// first we do a scan of the property
	openfpm::vector_gpu<aggregate<unsigned int>> idx;
	idx.resize(vd.size_local());

	auto ite = idx.getGPUIterator();

	create_index<<<ite.wthr,ite.thr>>>(idx.toKernel());

	// sort particles, so the particles to remove stay at the end
	mergesort((remove_type *)vd.getPropVector().template getDeviceBuffer<prp>(),(unsigned int *)idx.template getDeviceBuffer<0>(), idx.size(), mgpu::template less_t<remove_type>(), vd.getVC().getmgpuContext());

	openfpm::vector_gpu<aggregate<int>> mark;
	mark.resize(1);

	CudaMemory mem;
	mem.allocate(sizeof(int));
	mem.fill(0);

	// mark point, particle that stay and to remove
	find_buffer_offsets_no_prc<prp,decltype(vd.getPropVector().toKernel()),decltype(mark.toKernel())><<<ite.wthr,ite.thr>>>
			           (vd.getPropVector().toKernel(),(int *)mem.getDevicePointer(),mark.toKernel());

	mem.deviceToHost();

	// we have no particles to remove
	if (*(int *)mem.getPointer() == 0)
	{return;}

	// Get the mark point
	mark.template deviceToHost<0>();

	// than create an equivalent buffer prop and pos

	typename std::remove_reference<decltype(vd.getPosVector())>::type vd_pos_new;
	typename std::remove_reference<decltype(vd.getPropVector())>::type vd_prp_new;

	// resize them

	vd_pos_new.resize(mark.template get<0>(0));
	vd_prp_new.resize(mark.template get<0>(0));

	auto & vd_pos_old = vd.getPosVector();
	auto & vd_prp_old = vd.getPropVector();

	// now we copy from the old vector to the new one

	ite = vd_pos_old.getGPUIterator();

	copy_new_to_old<vector_type::dims><<<ite.wthr,ite.thr>>>(vd_pos_new.toKernel(),vd_prp_new.toKernel(),vd_pos_old.toKernel(),vd_prp_old.toKernel(),idx.toKernel());

	// and we swap

	vd.set_g_m(vd_pos_new.size());

	vd.getPosVector().swap(vd_pos_new);
	vd.getPropVector().swap(vd_prp_new);
}

#endif /* VECTOR_DIST_CUDA_FUNCS_CUH_ */
