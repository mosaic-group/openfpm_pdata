/*
 * ie_ghost_gpu.cuh
 *
 *  Created on: Aug 24, 2018
 *      Author: i-bird
 */

#ifndef IE_GHOST_GPU_CUH_
#define IE_GHOST_GPU_CUH_

#include "data_type/aggregate.hpp"

constexpr unsigned int lc_proc_ = 0;
constexpr unsigned int proc_ = 1;
constexpr unsigned int shift_id_ = 2;

template<unsigned int dim, typename T, typename cell_list_type, typename vb_int_box_type, typename vb_int_type>
__device__ __host__ inline unsigned int ghost_processorID_N_impl(const Point<dim,T> & p, cell_list_type & geo_cell, vb_int_box_type & vb_int_box, vb_int_type & vb_int)
{
	unsigned int cell = geo_cell.getCell(p);
	unsigned int sz = geo_cell.getNelements(cell);

	unsigned int n = 0;

	for (int i = 0 ; i < sz ; i++)
	{
		unsigned int bid = geo_cell.get(cell,i);

		if (Box<dim,T>(vb_int_box.get(bid)).isInsideNP(p) == true)
		{n++;}
	}

	return n;
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

	//! internal ghost box
	openfpm::vector_gpu_ker<aggregate<unsigned int,unsigned int,unsigned int>,layout_base> vb_int;

	// maximum rank
	unsigned int max_rank;

public:


	ie_ghost_gpu(CellList_cpu_ker<dim,T,Mem_fast_ker<Memory,memory_traits_lin,int>,shift<dim,T>> geo_cell,
			     openfpm::vector_gpu_ker<Box<dim, T>,layout_base> vb_int_box,
			     openfpm::vector_gpu_ker<aggregate<unsigned int,unsigned int,unsigned int>,layout_base> vb_int,
			     unsigned int max_rank)
	:geo_cell(geo_cell),vb_int_box(vb_int_box),vb_int(vb_int),max_rank(max_rank)
	{

	}

	ie_ghost_gpu(const ie_ghost_gpu<dim,T,Memory,layout_base> & ieg)
	:geo_cell(ieg.geo_cell),vb_int_box(ieg.vb_int_box),vb_int(ieg.vb_int),max_rank(ieg.max_rank)
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
		unsigned int cell = geo_cell.getCell(p);
		unsigned int sz = geo_cell.getNelements(cell);

		unsigned int n = 0;

		for (int i = 0 ; i < sz ; i++)
		{
			unsigned int bid = geo_cell.get(cell,i);

			if (Box<dim,T>(vb_int_box.get(bid)).isInsideNP(p) == true)
			{
				output.template get<0>(base+n) = vb_int.template get<proc_>(bid);
				output.template get<1>(base+n) = pi;

				n++;
			}
		}
	}

};



#endif /* IE_GHOST_GPU_CUH_ */
