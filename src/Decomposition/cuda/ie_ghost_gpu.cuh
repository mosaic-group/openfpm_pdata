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

template<unsigned int dim, typename T, typename cell_list_type, typename vb_int_box_type>
__device__ __host__ inline unsigned int ghost_processorID_N_impl(const Point<dim,T> & p, cell_list_type & geo_cell, vb_int_box_type & vb_int_proc)
{
	unsigned int cell = geo_cell.getCell(p);
	unsigned int sz = geo_cell.getNelements(cell);

	unsigned int n = 0;

	for (int i = 0 ; i < sz ; i++)
	{
		unsigned int bid = geo_cell.get(cell,i);

		unsigned int sz2 = vb_int_proc.template get<0>(bid).size();

		for (int j = 0 ; j < sz2 ; j++)
		{
			if (Box<dim,T>(vb_int_proc.template get<0>(bid).get(j)).isInsideNP(p) == true)
			{n++;}
		}
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
	openfpm::vector_gpu_ker<aggregate<openfpm::vector_gpu_ker<Box<dim, T>,layout_base>,int>,layout_base> vb_int_proc;

public:


	ie_ghost_gpu(CellList_cpu_ker<dim,T,Mem_fast_ker<Memory,memory_traits_lin,int>,shift<dim,T>> geo_cell,
				 openfpm::vector_gpu_ker<aggregate<openfpm::vector_gpu_ker<Box<dim, T>,layout_base>,int>,layout_base> vb_int_proc)
	:geo_cell(geo_cell),vb_int_proc(vb_int_proc)
	{

	}

	ie_ghost_gpu(const ie_ghost_gpu<dim,T,Memory,layout_base> & ieg)
	:geo_cell(ieg.geo_cell),vb_int_proc(ieg.vb_int_proc)
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
		return ghost_processorID_N_impl(p,geo_cell,vb_int_proc);
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

			unsigned int sz2 = vb_int_proc.template get<0>(bid).size();

			for (int j = 0 ; j < sz2 ; j++)
			{
				if (Box<dim,T>(vb_int_proc.template get<0>(bid).get(j)).isInsideNP(p) == true)
				{
					output.template get<0>(base+n) = vb_int_proc.template get<1>(bid);
					output.template get<1>(base+n) = pi;

					n++;
				}
			}
		}
	}

};



#endif /* IE_GHOST_GPU_CUH_ */
