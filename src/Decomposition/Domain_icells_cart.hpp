/*
 * Domain_icells_cart.hpp
 *
 *  Created on: Apr 27, 2019
 *      Author: i-bird
 */

#ifndef DOMAIN_ICELLS_CART_HPP_
#define DOMAIN_ICELLS_CART_HPP_

#include "Vector/map_vector.hpp"
#include "Space/Ghost.hpp"
#include "NN/CellList/CellList.hpp"
#include "NN/CellList/cuda/CellDecomposer_gpu_ker.cuh"
#include "Vector/map_vector_sparse.hpp"
#include <iomanip>

#ifdef __NVCC__

template<unsigned int dim, typename vector_sparse_type, typename CellDecomposer_type>
__global__ void insert_icell(vector_sparse_type vs, CellDecomposer_type cld, grid_key_dx<dim,int> start,grid_key_dx<dim,int> stop)
{
	vs.init();

	auto gk = grid_p<dim>::get_grid_point(cld.get_div_c());

	unsigned int b = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;

	bool out = false;
	for (unsigned int i = 0 ; i < dim ; i++)
	{
		gk.set_d(i,gk.get(i) + start.get(i));
		if (gk.get(i) > stop.get(i))
		{out = true;}
	}

	if (out == false)
	{
		auto id = cld.LinId(gk);

		vs.insert_b(id,b);
	}

	vs.flush_block_insert(b, threadIdx.x == 0 & threadIdx.y == 0 & threadIdx.z == 0 );
}

template<unsigned int dim, typename vector_sparse_type, typename CellDecomposer_type>
__global__ void insert_remove_icell(vector_sparse_type vs, vector_sparse_type vsi, CellDecomposer_type cld, grid_key_dx<dim,int> start,grid_key_dx<dim,int> stop)
{
	vs.init();
	vsi.init();

	auto gk = grid_p<dim>::get_grid_point(cld.get_div_c());

	unsigned int b = blockIdx.x + blockIdx.y * gridDim.x + blockIdx.z * gridDim.x * gridDim.y;

	bool out = false;
	for (unsigned int i = 0 ; i < dim ; i++)
	{
		gk.set_d(i,gk.get(i) + start.get(i));
		if (gk.get(i) > stop.get(i))
		{out = true;}
	}

	if (out == false)
	{
		auto id = cld.LinId(gk);

		vs.insert_b(id,b);
		vsi.remove_b(id,b);
	}

	vs.flush_block_insert(b, threadIdx.x == 0 & threadIdx.y == 0 & threadIdx.z == 0 );
	vsi.flush_block_remove(b, threadIdx.x == 0 & threadIdx.y == 0 & threadIdx.z == 0);
}

template<unsigned int dim, typename T, template<typename> class layout_base , typename Memory, typename cnt_type, typename ids_type, bool is_gpu>
struct CalculateInternalCells_impl
{
	template<typename VCluster_type>
	static void CalculateInternalCells(VCluster_type & v_cl,
			openfpm::vector<Box<dim,T>,Memory,layout_base> & ig_box,
			openfpm::vector<SpaceBox<dim,T>,Memory,layout_base> & domain,
			Box<dim,T> & pbox,
			T r_cut,
			const Ghost<dim,T> & enlarge,
			CellDecomposer_sm<dim,T,shift<dim,T>> & cd,
			openfpm::vector<aggregate<ids_type>,Memory,layout_base> & icells,
			openfpm::vector<aggregate<ids_type>,Memory,layout_base> & dcells)
	{

	}
};

template<unsigned int dim, typename T, template<typename> class layout_base , typename Memory, typename cnt_type, typename ids_type>
struct CalculateInternalCells_impl<dim,T,layout_base,Memory,cnt_type,ids_type,true>
{
	template<typename VCluster_type>
	static void CalculateInternalCells(VCluster_type & v_cl,
			openfpm::vector<Box<dim,T>,Memory,layout_base> & ig_box,
			openfpm::vector<SpaceBox<dim,T>,Memory,layout_base> & domain,
			Box<dim,T> & pbox,
			T r_cut,
			const Ghost<dim,T> & enlarge,
			CellDecomposer_sm<dim,T,shift<dim,T>> & cd,
			openfpm::vector<aggregate<ids_type>,Memory,layout_base> & icells,
			openfpm::vector<aggregate<ids_type>,Memory,layout_base> & dcells)
	{
#if 0

		// Division array
		size_t div[dim];

		// Calculate the parameters of the cell-list

		cl_param_calculate(pbox, div, r_cut, enlarge);

		openfpm::array<T,dim,cnt_type> spacing_c;
		openfpm::array<ids_type,dim,cnt_type> div_c;
		openfpm::array<ids_type,dim,cnt_type> off;

		for (size_t i = 0 ; i < dim ; i++)
		{
			spacing_c[i] = (pbox.getHigh(i) - pbox.getLow(i)) / div[i];
			off[i] = 1;
			// div_c must include offset
			div_c[i] = div[i] + 2*off[i];

		}

		shift_only<dim,T> t(Matrix<dim,T>::identity(),pbox.getP1());

		CellDecomposer_gpu_ker<dim,T,cnt_type,ids_type,shift_only<dim,T>> cld(spacing_c,div_c,off,t);
		grid_sm<dim,void> g = cld.getGrid();
		cd.setDimensions(pbox,div,off[0]);

		openfpm::vector_sparse_gpu<aggregate<unsigned int>> vs;
		openfpm::vector_sparse_gpu<aggregate<unsigned int>> vsi;

		vs.template setBackground<0>(0);

		// insert Domain cells

		for (size_t i = 0 ; i < domain.size() ; i++)
		{
			Box<dim,T> bx = SpaceBox<dim,T>(domain.get(i));

			auto pp2 = bx.getP2();

			for (size_t j = 0 ; j < dim ; j++)
			{pp2.get(j) = std::nextafter(pp2.get(j),pp2.get(j) - static_cast<T>(1.0));}

			auto p1 = cld.getCell(bx.getP1());
			auto p2 = cld.getCell(pp2);


			auto ite = g.getGPUIterator(p1,p2,256);

			if (ite.wthr.x == 0)
			{continue;}

			vsi.setGPUInsertBuffer(ite.nblocks(),256);

			CUDA_LAUNCH((insert_icell<dim>),ite,vsi.toKernel(),cld,ite.start,p2);

			vsi.template flush<>(v_cl.getgpuContext(),flush_type::FLUSH_ON_DEVICE);
		}

		// calculate the number of kernel launch

		for (size_t i = 0 ; i < ig_box.size() ; i++)
		{
			Box<dim,T> bx = ig_box.get(i);

			auto pp2 = bx.getP2();

			for (size_t j = 0 ; j < dim ; j++)
			{pp2.get(j) = std::nextafter(pp2.get(j),pp2.get(j) - static_cast<T>(1.0));}

			auto p1 = cld.getCell(bx.getP1());
			auto p2 = cld.getCell(pp2);

			auto ite = g.getGPUIterator(p1,p2,256);

			if (ite.wthr.x == 0)
			{continue;}

			vs.setGPUInsertBuffer(ite.nblocks(),256);
			vsi.setGPURemoveBuffer(ite.nblocks(),256);

			CUDA_LAUNCH(insert_remove_icell<dim>,ite,vs.toKernel(),vsi.toKernel(),cld,ite.start,p2);

			vs.template flush<>(v_cl.getgpuContext(),flush_type::FLUSH_ON_DEVICE);
			vsi.flush_remove(v_cl.getgpuContext(),flush_type::FLUSH_ON_DEVICE);
		}


		vs.swapIndexVector(icells);
		vsi.swapIndexVector(dcells);

#endif
	}
};

#endif

template<unsigned int dim, typename T, template<typename> class layout_base , typename Memory>
class domain_icell_calculator
{
	typedef unsigned int cnt_type;

	typedef int ids_type;

	openfpm::vector<aggregate<ids_type>,Memory,layout_base> icells;
	openfpm::vector<aggregate<ids_type>,Memory,layout_base> dcells;

	CellDecomposer_sm<dim,T,shift<dim,T>> cd;

	public:

	/*! \brief Calculate the subdomain that are in the skin part of the domain
	 *
       \verbatim

		+---+---+---+---+---+---+
		| 1 | 2 | 3 | 4 | 5 | 6 |
		+---+---+---+---+---+---+
		|28 |               | 7 |
		+---+               +---+
		|27 |               | 8 |
		+---+               +---+
		|26 |               | 9 |
		+---+   DOM1        +---+
		|25 |               |10 |
		+---+               +---+
		|24 |               |11 |
		+---+           +---+---+
		|23 |           |13 |12 |
		+---+-----------+---+---+
		|22 |           |14 |
		+---+           +---+
		|21 |   DOM2    |15 |
		+---+---+---+---+---+
		|20 |19 |18 | 17|16 |
		+---+---+---+---+---+    <----- Domain end here
                                      |
                        ^             |
                        |_____________|


       \endverbatim
	 *
	 * It does it on GPU or CPU
	 *
	 */
	template<typename VCluster_type>
	void CalculateInternalCells(VCluster_type & v_cl,
								openfpm::vector<Box<dim,T>,Memory,layout_base> & ig_box,
								openfpm::vector<SpaceBox<dim,T>,Memory,layout_base> & domain,
								Box<dim,T> & pbox,
								T r_cut,
								const Ghost<dim,T> & enlarge)
	{
#ifdef __NVCC__
		CalculateInternalCells_impl<dim,T,layout_base,Memory,cnt_type,ids_type,std::is_same<Memory,CudaMemory>::value>::CalculateInternalCells(v_cl,ig_box,domain,pbox,r_cut,enlarge,cd,icells,dcells);
#endif
	}

	/*! \brief Return the list of the internal cells
	 *
	 * \return the list of the internal cells
	 *
	 */
	openfpm::vector<aggregate<ids_type>,Memory,layout_base> & getIcells()
	{
		return icells;
	}

	/*! \brief Return the list of the internal cells
	 *
	 * \return the list of the internal cells
	 *
	 */
	openfpm::vector<aggregate<ids_type>,Memory,layout_base> & getDcells()
	{
		return dcells;
	}

	/*! \brief Given a cell index return the cell box
	 *
	 * \param ci cell index
	 *
	 * \return The box reppresenting the cell
	 */
	Box<dim,T> getBoxCell(unsigned int ci)
	{
		Box<dim,T> b;

		for (size_t i = 0 ; i < dim ; i++)
		{
			auto key = cd.getGrid().InvLinId(ci);
			Point<dim,T> p1 = cd.getOrig().get(i) - cd.getPadding(i)*cd.getCellBox().getHigh(i) ;

			b.setLow(i,p1.get(i) + key.get(i)*cd.getCellBox().getHigh(i));
			b.setHigh(i,p1.get(i) + ((key.get(i) + 1)*cd.getCellBox().getHigh(i)));
		}

		return b;
	}

	/*! \brief Get the grid base information about this cell decomposition
	 *
	 *
	 * \return the grid
	 *
	 */
	const grid_sm<dim,void> & getGrid()
	{
		return cd.getGrid();
	}
};


#endif /* DOMAIN_ICELLS_CART_HPP_ */
