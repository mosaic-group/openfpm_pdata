/*
 * space_key.hpp
 *
 *  Created on: Mar 6, 2015
 *      Author: Pietro Incardona
 */

#ifndef SPACE_KEY_HPP_
#define SPACE_KEY_HPP_


/*! \brief grid_key_dx is the key to access any element in the grid
 *
 * grid_key_dx is the key to access any element in the grid
 *
 * \tparam dim dimensionality of the grid
 * \tparam T type of the space float, double, complex, ...
 *
 */

template<unsigned int dim, typename T>
class space_key_dx
{
public:

	//! Constructor
	space_key_dx()
	{}

	//! Constructor from buffer
	space_key_dx(T (& k)[dim])
	{
		for (int i = 0 ; i < dim ; i++)
			this->k[i] = k[i];
	}

	//! Constructor from buffer
	space_key_dx(const T (& k)[dim])
	{
		for (int i = 0 ; i < dim ; i++)
			this->k[i] = k[i];
	}

	//! Construct a grid key from a list of numbers
	space_key_dx(const T v,const T...t)
	{
		k[dim-1] = v;
		invert_assign(t...);
	}

	/*! \brief Set to zero the key
	 */
	inline void zero()
	{
		for (int i = 0 ; i < dim ; i++)
		{
			k[i] = 0;
		}
	}

	/*! \brief Set to invalid the key
	 */
	inline void invalid()
	{
		for (int i = 0 ; i < dim ; i++)
		{
			k[i] = -1;
		}
	}

	/* \brief Check if two key are the same
	 *
	 * \param key_t key to check
	 *
	 * \return true if the two key are identical
	 *
	 */

	template<unsigned int dim_t> bool operator==(space_key_dx<dim_t,T> & key_t)
	{
		if (dim != dim_t)
		{
			return false;
		}

		// Check the two key index by index

		for (int i = 0 ; i < dim ; i++)
		{
			if (k[i] != key_t.k[i])
			{
				return false;
			}
		}

		// identical key
		return true;
	}

	//! set the grid key from a list of numbers
	template<typename a, typename ...T>void set(a v, T...t)
	{
		k[dim-1] = v;
		invert_assign(t...);
	}

	/*! \brief Get the i index
	 *
	 * \param i index to get
	 *
	 * \return the index value
	 *
	 */
	mem_id value(size_t i) const
	{
		return k[i];
	}

	/*! \brief Get the i index
	 *
	 *
	 * \param i index to get
	 *
	 * \return the index value
	 *
	 */
	mem_id get(size_t i) const
	{
		return k[i];
	}

	/*! \brief Set the i index
	 *
	 * Set the i index
	 *
	 * \param i index to set
	 * \param id value to set
	 *
	 */
	void set_d(size_t i, mem_id id)
	{
#ifdef DEBUG

		if (i >= dim)
			std::cerr << "grid_key_dx error: " << __FILE__ << " " << __LINE__ << " try to access dimension " << i << " on a grid_key_dx of size " << dim << "\n";

#endif
		k[i] = id;
	}

	//! structure that store all the index
	mem_id k[dim];

private:

	/*! \brief Recursively invert the assignment
	 *
	 * Recursively invert the assignment at compile-time
	 *
	 */
	template<typename a, typename ...T>void invert_assign(a v,T...t)
	{
		k[sizeof...(T)] = v;
		invert_assign(t...);
	}

	template<typename a, typename ...T>void invert_assign(a v)
	{
		k[0] = v;
	}

	void invert_assign()
	{
	}

};


#endif /* SPACE_KEY_HPP_ */
