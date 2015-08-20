/*
 * staggered_util.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_DIST_GRID_UTIL_HPP_
#define SRC_GRID_STAGGERED_DIST_GRID_UTIL_HPP_

/*! \brief Classes to get the number of components of the properties
 *
 */
template<typename T>
struct extends
{
	static inline size_t mul()
	{
		return 1;
	}

	static inline size_t dim()
	{
		return 0;
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1>
struct extends<T[N1]>
{
	static inline size_t mul()
	{
		return N1;
	}

	static inline size_t dim()
	{
		return 1;
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct extends<T[N1][N2]>
{
	static inline size_t mul()
	{
		return N1 * N2;
	}

	static inline size_t dim()
	{
		return 2;
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct extends<T[N1][N2][N3]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3;
	}

	static inline size_t dim()
	{
		return 3;
	}
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4>
struct extends<T[N1][N2][N3][N4]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4;
	}

	static inline size_t dim()
	{
		return 4;
	}
};

//! Partial specialization for N=5
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5>
struct extends<T[N1][N2][N3][N4][N5]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5;
	}

	static inline size_t dim()
	{
		return 5;
	}
};

//! Partial specialization for N=6
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6>
struct extends<T[N1][N2][N3][N4][N5][N6]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6;
	}

	static inline size_t dim()
	{
		return 6;
	}
};

//! Partial specialization for N=7
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7>
struct extends<T[N1][N2][N3][N4][N5][N6][N7]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7;
	}

	static inline size_t dim()
	{
		return 7;
	}
};

//! Partial specialization for N=8
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8>
struct extends<T[N1][N2][N3][N4][N5][N6][N7][N8]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8;
	}

	static inline size_t dim()
	{
		return 8;
	}
};

//! Partial specialization for N=9
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9>
struct extends<T[N1][N2][N3][N4][N5][N6][N7][N8][N9]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8 * N9;
	}

	static inline size_t dim()
	{
		return 9;
	}
};

//! Partial specialization for N=10
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9, size_t N10>
struct extends<T[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10]>
{
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8 * N9 * N10;
	}

	static inline size_t dim()
	{
		return 10;
	}
};

///////////////////// Staggered default positioning ////////////////////////

/*! \brief this class is a functor for "for_each" algorithm
 *
 * For each element of the boost::vector the operator() is called.
 * Is mainly used to produce a default position vector for each
 * property
 *
 * \tparam vector of properties
 *
 */

template<unsigned int dim, typename v, bool has_posMask>
struct stag_set_position
{
	openfpm::vector<comb<dim>> (& pos_prp)[boost::fusion::result_of::size<v>::type::value];

	//! It call the copy function for each property
	template<typename T>
	void operator()(T& t) const
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at<v,typename boost::mpl::int_<T::value>>::type prop;

		bool val = prop::stag_pos_mask[T::value];

		if (val == false)
			return;

		// Dimension of the object
		size_t dim_prp = extends<prop>::dim();

		// It is a scalar
		if (dim_prp == 0)
		{
			comb<dim> c;
			c.zero();

			// It stay in the center
			pos_prp[T::value].add(c);
		}
		else if (dim_prp == 1)
		{
			// It stay on the object of dimension dim-1 (Negative part)
			for (size_t i = 0 ; i < dim ; i++)
			{
				comb<dim> c;
				c.zero();
				c.value(i) = 1;

				pos_prp[T::value].add(c);
			}
		}
		else if (dim_prp == 2)
		{
			// Create an hypercube object
			HyperCube<dim> hyp;



			// Diagonal part live in
			for (size_t i = 0 ; i < dim ; i++)
			{
				comb<dim> c;
				c.zero();
				c.value(i) = 1;

				pos_prp[T::value].add(c);
			}
		}


	}
};


#endif /* SRC_GRID_STAGGERED_DIST_GRID_UTIL_HPP_ */
