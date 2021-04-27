/*
 * vector_dist_gpu.hpp
 *
 *  Created on: Jul 28, 2018
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_GPU_HPP_
#define VECTOR_DIST_GPU_HPP_

constexpr unsigned int POS_PROP = (unsigned int)-1;

#ifdef CUDA_GPU

#define GET_PARTICLE(vd) blockDim.x*blockIdx.x + threadIdx.x; if (blockDim.x*blockIdx.x + threadIdx.x >= vd.size_local()) {return;};
#define GET_PARTICLE_SORT(p,NN) if (blockDim.x*blockIdx.x + threadIdx.x >= NN.get_g_m()) {return;}\
							  else{p = NN.getDomainSortIds().template get<0>(blockDim.x*blockIdx.x + threadIdx.x);}


#define GET_PARTICLE_BY_ID(p,ids) if (blockDim.x*blockIdx.x + threadIdx.x >= ids.size()) {return;}\
							  else{p = ids.template get<0>(blockDim.x*blockIdx.x + threadIdx.x);}

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 *
 */
template<typename vector_dist_ker>
struct check_vector_dist_kernels
{
	//! op1
	const vector_dist_ker & o1;
	//! op2
	const vector_dist_ker & o2;

	bool check;

	/*! \brief constructor
	 *
	 * \param src source encapsulated object
	 * \param dst source encapsulated object
	 *
	 */
	inline check_vector_dist_kernels(const vector_dist_ker & o1, const vector_dist_ker & o2)
	:o1(o1),o2(o2),check(false)
	{};

	//! It call the copy function for each property
	template<typename T>
	__device__ __host__ inline void operator()(T& t)
	{
		check &= o1.template getPointer<T::value>() == o2.template getPointer<T::value>();
	}
};

template<unsigned int dim,
         typename St,
         typename prop,
         template<typename> class layout_base = memory_traits_inte>
class vector_dist_ker
{
	//! Ghost marker, all the particle with id > g_m are ghost all with g_m < are real particle
	int g_m = 0;

	//! Particle position vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	mutable openfpm::vector_gpu_ker<Point<dim,St>,layout_base> v_pos;

	//! Particle properties vector, (It has 2 elements) the first has real particles assigned to a processor
	//! the second element contain unassigned particles
	mutable openfpm::vector_gpu_ker<typename apply_transform<layout_base,prop>::type,layout_base> v_prp;

public:

	//! space type
	typedef St stype;

	//! dimensions of space
	static const unsigned int dims = dim;

	//! tag the type as a vector that run on kernel
	typedef int vector_kernel;

	//! Indicate this structure has a function to check the device pointer
	typedef int yes_has_check_device_pointer;

	vector_dist_ker()
	:g_m(0)
	{}

	vector_dist_ker(int g_m, const openfpm::vector_gpu_ker<Point<dim,St>,layout_base> & v_pos,
							 const openfpm::vector_gpu_ker<typename apply_transform<layout_base,prop>::type,layout_base> & v_prp)
	:g_m(g_m),v_pos(v_pos),v_prp(v_prp)
	{}

	/*! \brief return the number of particles (excluding ghost)
	 *
	 * \return the number of particles
	 *
	 */
	__device__ __host__ int size_local() const {return g_m;}

	/*! \brief return the number of particles
	 *
	 * \return the number of particles
	 *
	 */
	__device__ __host__ int size() const {return v_pos.size();}

	/*! \brief Get the position of an element
	 *
	 * see the vector_dist iterator usage to get an element key
	 *
	 * \param vec_key element
	 *
	 * \return the position of the element in space
	 *
	 */
	__device__ __host__ inline auto getPos(int vec_key) -> decltype(v_pos.template get<0>(vec_key))
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
	__device__ __host__ inline auto getPos(const vect_dist_key_dx & vec_key) -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
		return v_pos.template get<0>(vec_key.getKey());
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
	__device__ __host__ inline auto getPos(int vec_key) const -> decltype(v_pos.template get<0>(vec_key))
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
	__device__ __host__ inline auto getPos(const vect_dist_key_dx & vec_key) const -> decltype(v_pos.template get<0>(vec_key.getKey()))
	{
		return v_pos.template get<0>(vec_key.getKey());
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
	template<unsigned int id> __device__ __host__ inline auto getProp(int vec_key) -> decltype(v_prp.template get<id>(vec_key))
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
	template<unsigned int id> __device__ __host__ inline auto getProp(const vect_dist_key_dx & vec_key) -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
		return v_prp.template get<id>(vec_key.getKey());
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
	template<unsigned int id> __device__  __host__ inline auto getProp(int vec_key) const -> decltype(v_prp.template get<id>(vec_key))
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
	template<unsigned int id> __device__ __host__  inline auto getProp(const vect_dist_key_dx & vec_key) const -> decltype(v_prp.template get<id>(vec_key.getKey()))
	{
		return v_prp.template get<id>(vec_key.getKey());
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

	__host__ vector_dist_iterator getDomainIterator() const
	{
		std::cout << __FILE__ << ":" << __LINE__ << " error getDomainIterator used on a vector_dist_ker object is not allowed" << std::endl;

		return vector_dist_iterator(0,0);
	}

	/*! \brief Get an iterator that traverse the particles in the domain
	 *
	 * \return an iterator
	 *
	 */
	__host__ ite_gpu<1> getDomainIteratorGPU(size_t n_thr = default_kernel_wg_threads_) const
	{
		return v_pos.getGPUIteratorTo(g_m,n_thr);
	}

	/*! \brief Check that the two structures are the same (at level of pointers)
	 *
	 *
	 * \return
	 */
	__host__ bool operator==(const vector_dist_ker & v)
	{
		if (v.v_pos.template getPointer<0>() != v_pos.template getPointer<0>())
		{return false;}

		check_vector_dist_kernels<openfpm::vector_gpu_ker<prop,memory_traits_inte>> cv(this->v_prp,v.v_prp);

		cv.check = true;

		// Do the same for the properties
		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,prop::max_prop> >(cv);

		return cv.check;
	}

	/*! \brief This vector contain all particles so is not a subset
	 *
	 * \return false
	 *
	 */
	__host__ bool isSubset() const
	{
		return false;
	}

#ifdef SE_CLASS1

		/*! \brief Check if the device pointer is owned by this structure
		 *
		 * \return a structure pointer check with information about the match
		 *
		 */
		pointer_check check_device_pointer(void * ptr)
		{
			pointer_check pc;

			pc.match = false;

			// we check the position vector and the property vector
			pc = v_pos.check_device_pointer(ptr);

			if (pc.match == true)
			{
				pc.match_str = std::string("Particle index overflow in position (v_pos): ") + "\n" + pc.match_str;
				return pc;
			}

			pc = v_prp.check_device_pointer(ptr);
			if (pc.match == true)
			{
				pc.match_str = std::string("Particle index overflow in properties (v_prp): ") + "\n" + pc.match_str;
				return pc;
			}

			return pc;
		}

#endif
};

// This is a tranformation node for vector_distributed for the algorithm toKernel_tranform
template<template <typename> class layout_base, typename T>
struct toKernel_transform<layout_base,T,2>
{
	typedef typename apply_transform<layout_base,typename T::value_type>::type aggr;

	typedef vector_dist_ker<T::dims,typename T::stype,aggr,layout_base> type;
};

#endif

#endif /* VECTOR_DIST_GPU_HPP_ */
