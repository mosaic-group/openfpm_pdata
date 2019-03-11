/*
 * vector_dist_gpu.hpp
 *
 *  Created on: Jul 28, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_GPU_HPP_
#define VECTOR_DIST_GPU_HPP_

#ifdef CUDA_GPU

#define POS_PROP -1

#define GET_PARTICLE(vd) blockDim.x*blockIdx.x + threadIdx.x; if (blockDim.x*blockIdx.x + threadIdx.x >= vd.size()) {return;};
#define GET_PARTICLE_SORT(p,NN) if (blockDim.x*blockIdx.x + threadIdx.x >= NN.get_g_m()) {return;}\
							  else{p = NN.getDomainSortIds().template get<0>(blockDim.x*blockIdx.x + threadIdx.x);}

template<unsigned int dim,
         typename St,
         typename prop>
class vector_dist_ker
{
	//! Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	int g_m = 0;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	mutable openfpm::vector_gpu_ker<Point<dim,St>,memory_traits_inte> v_pos;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	mutable openfpm::vector_gpu_ker<prop,memory_traits_inte> v_prp;

public:

	//! space type
	typedef St stype;

	//! dimensions of space
	static const unsigned int dims = dim;

	vector_dist_ker(const openfpm::vector_gpu_ker<Point<dim,St>,memory_traits_inte> & v_pos, const openfpm::vector_gpu_ker<prop,memory_traits_inte> & v_prp)
	:v_pos(v_pos),v_prp(v_prp)
	{}

	/*! \brief return the number of particles (excluding ghost)
	 *
	 * \return the number of particles
	 *
	 */
	__device__ int size_local() {return g_m;}

	/*! \brief return the number of particles
	 *
	 * \return the number of particles
	 *
	 */
	__device__ int size() {return v_pos.size();}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	__device__ inline auto getPos(int vec_key) -> decltype(v_pos.template get<0>(vec_key))
	{
		return v_pos.template get<0>(vec_key);
	}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	__device__ inline auto getPos(int vec_key) const -> decltype(v_pos.template get<0>(vec_key))
	{
		return v_pos.template get<0>(vec_key);
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> __device__  inline auto getProp(int vec_key) -> decltype(v_prp.template get<id>(vec_key))
	{
		return v_prp.template get<id>(vec_key);
	}

	/*! \brief Get the property of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \tparam id property id
	 * \param vec_key vector element
	 *
	 * \return return the selected property of the vector element
	 *
	 */
	template<unsigned int id> __device__  inline auto getProp(int vec_key) const -> decltype(v_prp.template get<id>(vec_key))
	{
		return v_prp.template get<id>(vec_key);
	}

	/*! \brief Return the internal position vector
	 *
	 *
	 * \return the position vector
	 *
	 */
	__device__ openfpm::vector_gpu_ker<Point<dim,St>,memory_traits_inte> & getPosVector()
	{
		return v_pos;
	}

	/*! \return the internal property vector
	 *
	 * \return the property vector
	 *
	 */
	__device__ openfpm::vector_gpu_ker<prop,memory_traits_inte> & getPropVector()
	{
		return v_prp;
	}

};

// This is a tranformation node for vector_distributed for the algorithm toKernel_tranform
template<template <typename> class layout_base, typename T>
struct toKernel_transform<layout_base,T,2>
{
	typedef typename apply_transform<layout_base,typename T::value_type>::type aggr;

	typedef vector_dist_ker<T::dims,typename T::stype,aggr> type;
};

#endif

#endif /* VECTOR_DIST_GPU_HPP_ */
