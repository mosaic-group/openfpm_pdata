/*
 * CartDecomposition_gpu.hpp
 *
 *  Created on: Aug 7, 2018
 *      Author: i-bird
 */

#ifndef CARTDECOMPOSITION_GPU_HPP_
#define CARTDECOMPOSITION_GPU_HPP_

#include "ie_ghost_gpu.cuh"

template<unsigned int dim, typename bc_type, typename T2, typename fine_s_type, typename vsub_domain_type, typename box_type>
__device__ __host__ inline int processorID_impl(T2 & p,
												fine_s_type & fine_s,
												vsub_domain_type & sub_domains_global,
												const box_type & domain,
												const bc_type (& bc)[dim])
{
	// Get the number of elements in the cell

	int e = -1;
	int cl = fine_s.getCell(p);
	int n_ele = fine_s.getNelements(cl);

	int i = 0;
	for ( ; i < n_ele ; i++)
	{
		e = fine_s.get(cl,i);

		if (sub_domains_global.template get<0>(e).isInsideNP_with_border(p,domain,bc) == true)
		{
			break;
		}
	}

#if defined(SE_CLASS1)

	if (n_ele == 0)
	{
		printf("CartDecomposition_gpu.cuh:processorID_impl, error I cannot detect in which processor this particle go \n");
		return -1;
	}

	if (i == n_ele)
	{
		printf("CartDecomposition_gpu.cuh:processorID_impl, error I cannot detect in which processor this particle go because of round-off inconsistencies \n");
		return -1;
	}

#endif


	/* coverty[negative_returns] */
	return sub_domains_global.template get<1>(e);
}

/*! \brief Apply boundary condition to the point
 *
 * If the particle go out to the right, bring back the particle on the left
 * in case of periodic, nothing in case of non periodic
 *
 * \param pt encapsulated point object (it's coordinated are changed according the
 *        the explanation before)
 *
 */
template<unsigned int dim, typename St, typename Mem> inline __device__ void applyPointBC_no_dec(Box<dim,St> & domain, periodicity_int<dim> & bc, encapc<1,Point<dim,St>,Mem> && pt)
{
	for (size_t i = 0 ; i < dim ; i++)
	{
		if (bc.bc[i] == PERIODIC)
		{pt.template get<0>()[i] = openfpm::math::periodic_l(pt.template get<0>()[i],domain.getHigh(i),domain.getLow(i));}
	}
}

template<unsigned int dim, typename T, typename Memory, template <typename> class layout_base>
class CartDecomposition_gpu : public ie_ghost_gpu<dim,T,Memory,layout_base>
{
	CellList_cpu_ker<dim,T,Mem_fast_ker<Memory,memory_traits_lin,int>,shift<dim,T>> clk;

	Box<dim,T> domain;

	int bc[dim];

	openfpm::vector_gpu_ker<Box_map<dim, T>,layout_base> sub_domains_global;

	/*! \brief Apply boundary condition to the point
	 *
	 * If the particle go out to the right, bring back the particle on the left
	 * in case of periodic, nothing in case of non periodic
	 *
	 * \param pt Point to apply the boundary conditions.(it's coordinated are changed according the
	 *        the explanation before)
	 *
	 */
	__device__ __host__ void applyPointBC(Point<dim,T> & pt) const
	{
		for (int i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
			{pt.get(i) = openfpm::math::periodic_l(pt.get(i),domain.getHigh(i),domain.getLow(i));}
		}
	}

public:

	CartDecomposition_gpu(CellList_cpu_ker<dim,T,Mem_fast_ker<Memory,memory_traits_lin,int>,shift<dim,T>> clk,
						  ie_ghost_gpu<dim,T,Memory,layout_base> ieg,
						  openfpm::vector_gpu_ker<Box_map<dim, T>,layout_base> sub_domains_global,
						  const Box<dim,T> & domain,
						  const int (& bc)[dim])
	:ie_ghost_gpu<dim,T,Memory,layout_base>(ieg),clk(clk),domain(domain),sub_domains_global(sub_domains_global)
	{
		for (int s = 0 ; s < dim ; s++)
		{this->bc[s] = bc[s];}
	}

	CartDecomposition_gpu(const CartDecomposition_gpu<dim,T,Memory,layout_base> & dec)
	:ie_ghost_gpu<dim,T,Memory,layout_base>(dec),clk(dec.clk),domain(dec.domain),sub_domains_global(dec.sub_domains_global)
	{
		for (int s = 0 ; s < dim ; s++)
		{this->bc[s] = dec.bc[s];}
	}

	/*! \brief Given a point return in which processor the point/particle should go
	 *
	 * Boundary conditions are considered
	 *
	 * \param p point
	 *
	 * \return processorID
	 *
	 */
	__device__ __host__ int inline processorIDBC(const Point<dim,T> & p)
	{
		Point<dim,T> pt = p;
		this->applyPointBC(pt);

		return processorID_impl(pt,clk,sub_domains_global,domain,bc);
	}

	/*! \brief Apply boundary condition to the point
	 *
	 * If the particle go out to the right, bring back the particle on the left
	 * in case of periodic, nothing in case of non periodic
	 *
	 * \param pt encapsulated point object (it's coordinated are changed according the
	 *        the explanation before)
	 *
	 */
	template<typename Mem> __device__ __host__ void applyPointBC(encapc<1,Point<dim,T>,Mem> && pt) const
	{
		for (size_t i = 0 ; i < dim ; i++)
		{
			if (bc[i] == PERIODIC)
			{pt.template get<0>()[i] = openfpm::math::periodic_l(pt.template get<0>()[i],domain.getHigh(i),domain.getLow(i));}
		}
	}


	/*! \brief Given a point return in which processor the particle should go
	 *
	 * \param p point
	 *
	 * \return processorID
	 *
	 */
	__device__ __host__ int inline processorID(const Point<dim,T> &pt)
	{
		return processorID_impl(pt,clk,sub_domains_global,domain,bc);
	}
};

#endif /* CARTDECOMPOSITION_GPU_HPP_ */
