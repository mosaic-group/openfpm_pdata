/*
 * vector_dist_cuda_funcs.cuh
 *
 *  Created on: Aug 14, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_CUDA_FUNCS_CUH_
#define VECTOR_DIST_CUDA_FUNCS_CUH_

#include "Vector/util/vector_dist_funcs.hpp"
#include "Decomposition/common.hpp"
#include "lib/pdata.hpp"
#include "util/cuda/kernels.cuh"
#include "util/cuda/scan_ofp.cuh"
#include "util/cuda/reduce_ofp.cuh"
#include "memory/CudaMemory.cuh"

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

	vd_prp.template set<prp ...>(p,vd_prp_ord,nss.template get<0>(p));
}

template<typename vector_pos_type, typename vector_prp_type, typename stns_type, unsigned int ... prp>
__global__ void merge_sort_all(vector_pos_type vd_pos, vector_prp_type vd_prp,
		                        vector_pos_type v_pos_ord, vector_prp_type vd_prp_ord,
		                        stns_type nss)
{
	int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd_pos.size()) return;

	vd_pos.template set<0>(p,v_pos_ord,nss.template get<0>(p));

	vd_prp.set(p,vd_prp_ord,nss.template get<0>(p));
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


template<typename vector_prp_type_out, typename vector_prp_type_in, unsigned int ... prp>
__global__  void process_ghost_particles_prp_put(vector_prp_type_out m_prp,
		     	 	 	 	 	 	   	     vector_prp_type_in  v_prp, unsigned int offset)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= m_prp.size()) return;

    process_ghost_device_particle_prp<vector_prp_type_out,vector_prp_type_in,prp...>(i,offset,m_prp,v_prp);
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

template<unsigned int dim, typename St, typename vector_of_box, typename vector_of_shifts, typename vector_type,  typename output_type>
__global__ void num_shift_ghost_each_part(vector_of_box box_f, vector_of_shifts box_f_sv, vector_type vd,  output_type out, unsigned int g_m)
{
	unsigned int old_shift = (unsigned int)-1;
	int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= g_m) return;

    Point<dim,St> xp = vd.template get<0>(p);

    unsigned int n = 0;

    for (unsigned int i = 0 ; i < box_f.size() ; i++)
    {
    	unsigned int shift_actual = box_f_sv.template get<0>(i);
    	bool sw = (old_shift == shift_actual)?true:false;

    	if (Box<dim,St>(box_f.get(i)).isInsideNP(xp) == true && sw == false)
    	{
    		old_shift = shift_actual;
    		n++;
    	}
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
		                              output_type output, unsigned int offset,unsigned int g_m)
{
	unsigned int old_shift = (unsigned int)-1;
	int p = threadIdx.x + blockIdx.x * blockDim.x;

    if (p >= g_m) return;

    Point<dim,St> xp = v_pos.template get<0>(p);

    unsigned int base_o = start.template get<0>(p);
    unsigned int base = base_o + offset;


    unsigned int n = 0;

    for (unsigned int i = 0 ; i < box_f.size() ; i++)
    {
    	unsigned int shift_actual = box_f_sv.template get<0>(i);
    	bool sw = (old_shift == shift_actual)?true:false;

    	if (Box<dim,St>(box_f.get(i)).isInsideNP(xp) == true && sw == false)
    	{

#pragma unroll
    		for (unsigned int j = 0 ; j < dim ; j++)
    		{
    			v_pos.template get<0>(base+n)[j] = xp.get(j) - shifts.template get<0>(shift_actual)[j];
    		}

    		output.template get<0>(base_o+n) = p;
    		output.template get<1>(base_o+n) = shift_actual;

    		v_prp.set(base+n,v_prp.get(p));

    		old_shift = shift_actual;
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
struct _add_: gpu::plus_t<red_type>
{};

template<typename red_type>
struct _max_: gpu::maximum_t<red_type>
{};

template<unsigned int prp, template <typename> class op, typename vector_type>
auto reduce_local(vector_type & vd) -> typename std::remove_reference<decltype(vd.template getProp<prp>(0))>::type
{
	typedef typename std::remove_reference<decltype(vd.template getProp<prp>(0))>::type reduce_type;

	CudaMemory mem;
	mem.allocate(sizeof(reduce_type));

	openfpm::reduce((reduce_type *)vd.getPropVector(). template getDeviceBuffer<prp>(),
			            vd.size_local(), (reduce_type *)mem.getDevicePointer() ,
			            op<reduce_type>(), vd.getVC().getgpuContext());

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

template<unsigned int dim, typename vector_pos_type, typename vector_prp_type, typename ids_type>
__global__ void copy_new_to_old(vector_pos_type vd_pos_dst, vector_prp_type vd_prp_dst, vector_pos_type vd_pos_src, vector_prp_type vd_prp_src, ids_type idx)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= vd_prp_dst.size()) return;

    for (unsigned int k = 0 ; k < dim ; k++)
    {vd_pos_dst.template get<0>(i)[k] = vd_pos_src.template get<0>(idx.template get<0>(i))[k];}

    vd_prp_dst.set(i,vd_prp_src,idx.template get<0>(i));
}

template<unsigned int dim, unsigned int prp, typename vector_pos_type, typename vector_prp_type, typename scan_type>
__global__ void copy_new_to_old_by_scan(vector_pos_type vd_pos_dst, vector_prp_type vd_prp_dst, vector_pos_type vd_pos_src, vector_prp_type vd_prp_src, scan_type scan)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= scan.size()) return;

    auto sc = scan.template get<0>(i);

    if (vd_prp_src.template get<prp>(i) == 0) return;

    for (unsigned int k = 0 ; k < dim ; k++)
    {vd_pos_dst.template get<0>(sc)[k] = vd_pos_src.template get<0>(i)[k];}

    vd_prp_dst.set(sc,vd_prp_src,i);
}


template<unsigned int prp,typename vector_type>
__global__ void flip_one_to_zero(vector_type vd)
{
    int i = threadIdx.x + blockIdx.x * blockDim.x;

    if (i >= vd.size_local()) return;

    vd.template getProp<prp>(i) = (vd.template getProp<prp>(i) == 0);
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
void remove_marked(vector_type & vd, const int n = 1024)
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

	if (vd.size_local() == 0)
	{return;}

	typedef typename boost::mpl::at<typename vector_type::value_type::type,boost::mpl::int_<prp>>::type remove_type;

	// because we mark the one to remove we flip the one to zero and the zeros to one

	auto ite = vd.getDomainIteratorGPU(n);

	CUDA_LAUNCH((flip_one_to_zero<prp>),ite,vd.toKernel());

	// first we scan

	openfpm::vector_gpu<aggregate<remove_type>> idx;

	if (mem_tmp.ref() == 0)
	{mem_tmp.incRef();}

	idx.setMemory(mem_tmp);
	idx.resize(vd.size_local());

	openfpm::scan((remove_type *)vd.getPropVector().template getDeviceBuffer<prp>(),vd.size_local(),(remove_type *)idx.template getDeviceBuffer<0>(),vd.getVC().getgpuContext());

	// Check if we marked something

	idx.template deviceToHost<0>(idx.size()-1,idx.size()-1);
	vd.template deviceToHostProp<prp>(vd.size_local()-1,vd.size_local()-1);

	int n_marked = vd.size_local() - (vd.template getProp<prp>(vd.size_local()-1) + idx.template get<0>(idx.size()-1));

	if (n_marked == 0)
	{
		// No particle has been marked Nothing to do

		return;
	}

	/////////////////////////////////

	typename std::remove_reference<decltype(vd.getPosVector())>::type vd_pos_new;
	typename std::remove_reference<decltype(vd.getPropVector())>::type vd_prp_new;

	// we reuse memory. this give us the possibility to avoid allocation and make the remove faster

	// Reference counter must be set to zero

/*	for (int i = 0 ; i < MAX_NUMER_OF_PROPERTIES ; i++)
	{
		for (int j = 0 ; j < exp_tmp2[i].ref() ; j++)
		{exp_tmp2[i].decRef();}
	}

	for (int j = 0 ; j < exp_tmp.ref() ; j++)
	{exp_tmp.decRef();}*/

	//vd_pos_new.setMemory(exp_tmp);
	//vd_prp_new.setMemoryArray((CudaMemory *)&exp_tmp2);

	// resize them

	vd_pos_new.resize(vd.size_local() - n_marked);
	vd_prp_new.resize(vd.size_local() - n_marked);

	auto & vd_pos_old = vd.getPosVector();
	auto & vd_prp_old = vd.getPropVector();

	CUDA_LAUNCH((copy_new_to_old_by_scan<vector_type::dims,prp>),ite,vd_pos_new.toKernel(),vd_prp_new.toKernel(),vd_pos_old.toKernel(),vd_prp_old.toKernel(),idx.toKernel());

	vd.set_g_m(vd_pos_new.size());

	vd.getPosVector().swap_nomode(vd_pos_new);
	vd.getPropVector().swap_nomode(vd_prp_new);

	// Increment v_pos_new and vd_prp_new memory counters

	vd.setReferenceCounterToOne();
}

template<unsigned int prp, typename functor, typename particles_type, typename out_type>
__global__ void mark_indexes(particles_type vd, out_type out, unsigned int g_m)
{
	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= vd.size())	{return;}

	out.template get<0>(p) = functor::check(vd.template get<prp>(p)) == true && p < g_m;
}

template<typename out_type, typename ids_type>
__global__ void fill_indexes(out_type scan, ids_type ids)
{
	unsigned int p = threadIdx.x + blockIdx.x * blockDim.x;

	if (p >= scan.size()-1)	{return;}

	auto sp = scan.template get<0>(p);
	auto spp = scan.template get<0>(p+1);

	if (sp != spp)
	{ids.template get<0>(scan.template get<0>(p)) = p;}
}

/*! \brief get the particle index that satify the functor condition
 *
 * This function can be used to collect the indexes of the particles of a particular type.
 * Write a functor that return true when a particle of a particular type is identified
 * and ids will contain the indexes for which the functor return true.
 *
 * \tparam prp property to pass to the functor
 *
 * \param vd distributed vector
 *
 */
template<unsigned int prp, typename functor, typename vector_type, typename ids_type>
void get_indexes_by_type(vector_type & vd, ids_type & ids, size_t end ,gpu::ofp_context_t & context)
{
	// first we do a scan of the property
	openfpm::vector_gpu<aggregate<unsigned int>> scan;

	scan.setMemory(mem_tmp);
	scan.resize(vd.size()+1);

	auto ite = scan.getGPUIterator();

	CUDA_LAUNCH((mark_indexes<prp,functor>),ite,vd.toKernel(),scan.toKernel(),end);

	openfpm::scan((unsigned int *)scan.template getDeviceBuffer<0>(),scan.size(),(unsigned int *)scan.template getDeviceBuffer<0>(),context);

	// get the number of marked particles
	scan.template deviceToHost<0>(scan.size()-1,scan.size()-1);
	size_t nf = scan.template get<0>(scan.size()-1);
	ids.resize(nf);

	CUDA_LAUNCH(fill_indexes,ite,scan.toKernel(),ids.toKernel());
}

#endif /* VECTOR_DIST_CUDA_FUNCS_CUH_ */
