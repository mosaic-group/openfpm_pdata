/*
 * vector_dist_operators_list_ker.hpp
 *
 *  Created on: Jun 7, 2019
 *      Author: i-bird
 */

#ifndef VECTOR_DIST_OPERATORS_LIST_KER_HPP_
#define VECTOR_DIST_OPERATORS_LIST_KER_HPP_

template<typename T>
struct ref_wrap
{
	T & v;

	ref_wrap(T & v)
	:v(v) {}

	ref_wrap & operator=(const ref_wrap<T> & rw)
	{
		v = rw.v;
		return *this;
	}
};

/*! \brief This class contain a list of all tracked vector_dist_ker around.
 *
 * In short suppose to do a auto vdk = vd.toKernel() vdk wrap the cuda pointer of vd. On the other hand
 * if in vd we use instruction that produce reallocations vdk contain invalid pointers. this class contail a
 * list of all vector_dist_kernel, such that if a reallocation happen the pointer are updated on the vector_dist_kernel
 *
 */
template<typename vector_dist_ker_type>
class vector_dist_ker_list
{
	openfpm::vector<ref_wrap<vector_dist_ker_type>> vkers;

public:

	/*! \brief Add a new vector_dist_kernel to track
	 *
	 * \param v vector_dist_kernel to track
	 *
	 *
	 */
	void add(vector_dist_ker_type & v)
	{
		ref_wrap<vector_dist_ker_type> rw(v);

		vkers.add(rw);

		if (vkers.size() >= 64)
		{
			std::cout << __FILE__ << ":" << __LINE__ << " The array of tracked vector_dist_ker become suspiciously big, are there memory leak ? " << std::endl;
		}
	}

	/*! \brief Update the addresses of all vector_dist_kernels around
	 *
	 * \param v vector_dist_kernel to track
	 *
	 *
	 */
	void update(const vector_dist_ker_type & v)
	{
		for (size_t i = 0 ; i < vkers.size() ; i++)
		{
			vkers.get(i).v = v;
		}
	}

	/*! \brief Remove one vector_dist_kernels entry
	 *
	 * \param v vector_dist_kernel to remove
	 *
	 *
	 */
	void remove(vector_dist_ker_type & v)
	{
		for (size_t i = 0 ; i < vkers.size() ; i++)
		{
			if (&vkers.get(i).v == &v)
			{
				vkers.remove(i);
				break;
			}
		}
	}

	/*! \brief Return the number of entries
	 *
	 * \return the number of entries
	 *
	 */
	size_t n_entry()
	{
		return vkers.size();
	}

	/*! \brief Check that all the entries are aligned to the latest vector_dist_ker_type
	 *
	 * \return true if all entries are alligned
	 *
	 */
	bool check(const vector_dist_ker_type & v)
	{
		for (size_t i = 0 ; i < vkers.size() ; i++)
		{
			if (!(vkers.get(i).v == v))
			{
				return false;
			}
		}

		return true;
	}
};

#endif /* VECTOR_DIST_OPERATORS_LIST_KER_HPP_ */
