/*
 * ie_ghost_gpu.cuh
 *
 *  Created on: Aug 24, 2018
 *      Author: i-bird
 */

#ifndef IE_GHOST_GPU_CUH_
#define IE_GHOST_GPU_CUH_

#include "Decomposition/common.hpp"
#include "data_type/aggregate.hpp"
#include "Space/Shape/Point.hpp"
#include "Space/Shape/Box.hpp"
#include "NN/CellList/cuda/CellList_cpu_ker.cuh"
#include "NN/Mem_type/MemFast.hpp"
#include "NN/CellList/CellDecomposer.hpp"

constexpr unsigned int lc_proc_ = 0;
constexpr unsigned int proc_ = 1;
constexpr unsigned int shift_id_ = 2;

template<typename output_type>
struct ID_operation
{
	output_type & output;

	__device__ __host__ ID_operation(output_type & output)
	:output(output)
	{}

	__device__ __host__ inline void op(unsigned int base, unsigned int n, unsigned int proc_act, unsigned int shift_act, unsigned int pi)
	{
		output.template get<0>(base + n) = proc_act;
		output.template get<1>(base + n) = (unsigned long int)pi + (((unsigned long int)shift_act) << 32);
	}
};

struct N_operation
{
	__device__ __host__ inline void op(unsigned int base, unsigned int n, unsigned int proc_act, unsigned int shift_act, unsigned int pi)
	{
	}
};

template<unsigned int dim, typename T, typename cell_list_type, typename vb_int_box_type, typename vb_int_type, typename operation>
__device__ __host__ inline unsigned int ghost_processorID_general_impl(const Point<dim,T> & p,
																 unsigned int base,
																 unsigned int pi,
																 cell_list_type & geo_cell,
																 vb_int_box_type & vb_int_box,
																 vb_int_type & vb_int,
																 operation & op)
{
	unsigned int cell = geo_cell.getCell(p);
	unsigned int sz = geo_cell.getNelements(cell);

	unsigned int n = 0;

	bool switch_prc = false;

	if (sz != 0)
	{
		int i = 0;
		unsigned int bid = geo_cell.get(cell,0);
		unsigned int proc_prev = vb_int.template get<proc_>(bid);
		unsigned int shift_prev = vb_int.template get<shift_id_>(bid);
		unsigned int proc_act;
		unsigned int shift_act;

		if (Box<dim,T>(vb_int_box.get(bid)).isInsideNP(p) == true)
		{
			op.op(base,n,proc_prev,shift_prev,pi);

			switch_prc = true;
			n++;
		}

		i++;

		for ( ; i < sz ; i++)
		{
			unsigned int bid = geo_cell.get(cell,i);
			proc_act = vb_int.template get<proc_>(bid);
			shift_act = vb_int.template get<shift_id_>(bid);

			switch_prc = (proc_act == proc_prev && shift_act == shift_prev) & switch_prc;

			if (Box<dim,T>(vb_int_box.get(bid)).isInsideNP(p) == true && switch_prc == false)
			{
				op.op(base,n,proc_act,shift_act,pi);

				switch_prc = true;
				n++;
			}
			proc_prev = proc_act;
			shift_prev = shift_act;
		}
	}

	return n;
}

template<unsigned int dim, typename T, typename cell_list_type, typename vb_int_box_type, typename vb_int_type>
__device__ __host__ inline unsigned int ghost_processorID_N_impl(const Point<dim,T> & p,
																 cell_list_type & geo_cell,
																 vb_int_box_type & vb_int_box,
																 vb_int_type & vb_int)
{
	N_operation op;

	return ghost_processorID_general_impl(p,0,0,geo_cell,vb_int_box,vb_int,op);
}

/*! \brief structure that store and compute the internal and external local ghost box. Version usable in kernel
 *
 * \tparam dim is the dimensionality of the physical domain we are going to decompose.
 * \tparam T type of the space we decompose, Real, Integer, Complex ...
 *
 * \see CartDecomposition
 *
 */
template<unsigned int dim, typename T, typename Memory, template<typename> class layout_base>
class ie_ghost_gpu
{

	//! Cell-list that store the geometrical information of the internal ghost boxes
	CellList_cpu_ker<dim,T,Mem_fast_ker<Memory,memory_traits_lin,int>,shift<dim,T>> geo_cell;

	//! internal ghost box
	openfpm::vector_gpu_ker<Box<dim, T>,layout_base> vb_int_box;

	//! internal ghost box processor infos
	openfpm::vector_gpu_ker<aggregate<unsigned int,unsigned int,unsigned int>,layout_base> vb_int;

public:


	ie_ghost_gpu(CellList_cpu_ker<dim,T,Mem_fast_ker<Memory,memory_traits_lin,int>,shift<dim,T>> geo_cell,
				 openfpm::vector_gpu_ker<Box<dim, T>,layout_base> vb_int_box,
				 openfpm::vector_gpu_ker<aggregate<unsigned int,unsigned int,unsigned int>,layout_base> vb_int)
	:geo_cell(geo_cell),vb_int_box(vb_int_box),vb_int(vb_int)
	{

	}

	ie_ghost_gpu(const ie_ghost_gpu<dim,T,Memory,layout_base> & ieg)
	:geo_cell(ieg.geo_cell),vb_int_box(ieg.vb_int_box),vb_int(ieg.vb_int)
	{}

	/*! \brief Get the cell from the particle position
	 *
	 * \param p position of the particle
	 *
	 */
	__device__ inline unsigned int ghost_processorID_cell(const Point<dim,T> & p)
	{
		return geo_cell.getCell(p);
	}

	/*! \brief Get the number of processor a particle must sent
	 *
	 * \param p position of the particle
	 *
	 */
	__device__ inline unsigned int ghost_processorID_N(const Point<dim,T> & p)
	{
		return ghost_processorID_N_impl(p,geo_cell,vb_int_box,vb_int);
	}

	/*! \brief Get the number of processor a particle must sent
	 *
	 * \param p position of the particle
	 *
	 */
	template<typename output_type> __device__ inline void ghost_processor_ID(const Point<dim,T> & p, output_type & output, unsigned int base, unsigned int pi)
	{
		ID_operation<output_type> op(output);

		ghost_processorID_general_impl(p,base,pi,geo_cell,vb_int_box,vb_int,op);
	}

};



#endif /* IE_GHOST_GPU_CUH_ */
