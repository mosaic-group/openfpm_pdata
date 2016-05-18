/*
 * staggered_grid_dist_copy.hpp
 *
 *  Created on: Apr 9, 2016
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_DIST_GRID_COPY_HPP_
#define SRC_GRID_STAGGERED_DIST_GRID_COPY_HPP_

/*!	\brief Add scalar elements
 *
 * \tparam copy_type Type that should be copied
 * \tparam T property id to copy
 * \tparam Grid_src Staggered source Grid
 * \tparam sa dimensionality of the array 0 is a scalar
 *
 */
template<typename copy_type, typename Tsrc, typename Tdst, typename Grid_src, typename Grid_dst, int sa>
struct interp_ele_sca_array
{
	inline static void interp(Grid_dst & grid_dst, const grid_dist_key_dx<Grid_dst::dims> & key_dst ,const Grid_src & x, const grid_dist_key_dx<Grid_src::dims> & key_src, const openfpm::vector<std::vector<comb<Grid_src::dims>>> & interp_pos)
	{
		typedef typename boost::mpl::at<Tdst,Tsrc>::type Tdst_ele;

		copy_type division = 0.0;

		for (size_t i = 0 ; i < interp_pos.get(0).size() ; i++)
		{
				auto key_m = key_src.move(interp_pos.get(0)[i]);

				grid_dst.template get<Tdst_ele::value>(key_dst) += x.template get<Tsrc::value>(key_m);

				division += 1.0;
		}
		grid_dst.template get<Tdst_ele::value>(key_dst) /= division;
	}
};

/*! \brief Add 1D array elements
 *
 * spacialization in case of 1D array
 *
 * \tparam copy_type Type that should be copied
 * \tparam T property id to copy
 * \tparam Ev Type of source the Vector
 *
 */
template<typename copy_type, typename Tsrc, typename Tdst, typename Grid_src, typename Grid_dst>
struct interp_ele_sca_array<copy_type,Tsrc,Tdst,Grid_src,Grid_dst,1>
{
	inline static void interp(Grid_dst & grid_dst,
								const grid_dist_key_dx<Grid_dst::dims> & key_dst ,
								const Grid_src & x,
								const grid_dist_key_dx<Grid_src::dims> & key_src,
								const openfpm::vector<std::vector<comb<Grid_src::dims>>> & interp_pos)
	{
		typename std::remove_all_extents<copy_type>::type division;
		typedef typename boost::mpl::at<Tdst,Tsrc>::type Tdst_ele;

		for (size_t j = 0 ; j < std::extent<copy_type>::value ; j++)
		{
			division = 0.0;
			for (size_t i = 0 ; i < interp_pos.get(j).size() ; i++)
			{
				auto key_m = key_src.move(interp_pos.get(j)[i]);

				grid_dst.template get<Tdst_ele::value>(key_dst)[j] += x.template get<Tsrc::value>(key_m)[j];

				division += 1.0;
			}
			grid_dst.template get<Tsrc::value>(key_dst)[j] /= division;
		}
	}
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to interpolate from the staggered grid to the normal target grid
 *
 * \tparam Tdst destination property of the normal grid
 * \tparam Type of the destination grid
 * \tparam Type of the source grid
 *
 */
template<typename Tdst, typename Grid_dst, typename Grid_src, unsigned int nst_pos>
struct interp_ele
{
	//! destination point
	const grid_dist_key_dx<Grid_dst::dims> key_dst;

	//! destination grid
	Grid_dst & grid_dst;

	//! source point
	grid_dist_key_dx<Grid_dst::dims> key_src;

	//! For each properties [] for each components (openfpm::vector) interpolants points positions (std::vector<comb>)
	openfpm::vector<std::vector<comb<Grid_dst::dims>>> (&interp_pos)[nst_pos];

	//! source grid
	const Grid_src & x;

	/*! \brief constructor
	 *
	 * It define the interpolation parameters.
	 *
	 * \param key_dst destination point
	 * \param grid_dst Destination grid
	 * \param key_src source point
	 * \param grid_src Source grid
	 *
	 */
	inline interp_ele(const grid_dist_key_dx<Grid_dst::dims> & key_dst, Grid_dst & grid_dst, const Grid_src & x, const grid_dist_key_dx<Grid_src::dims> & key_src, openfpm::vector<std::vector<comb<Grid_src::dims>>> (&interp_pos)[nst_pos])
	:key_dst(key_dst),grid_dst(grid_dst),key_src(key_src),interp_pos(interp_pos),x(x){};


#ifdef SE_CLASS1
	/*! \brief Constructor
	 *
	 * Calling this constructor produce an error. This class store the reference of the object,
	 * this mean that the object passed must not be a temporal object
	 *
	 */
	inline interp_ele(const grid_dist_key_dx<Grid_dst::dims> & key_dst, Grid_dst && grid_dst, const Grid_src & x, const grid_dist_key_dx<Grid_src::dims> & key_src, openfpm::vector<std::vector<comb<Grid_src::dims>>> (&interp_pos)[nst_pos])
	:key_dst(key_dst),grid_dst(grid_dst),key_src(key_src),interp_pos(interp_pos),x(x)
	{std::cerr << "Error: " <<__FILE__ << ":" << __LINE__ << " Passing a temporal object";};
#endif

	//! It call the copy function for each property
	template<typename Tsrc>
	inline void operator()(Tsrc& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at_c<typename Grid_dst::value_type::type,Tsrc::value>::type copy_type;

		interp_ele_sca_array<copy_type,Tsrc,Tdst,Grid_src,Grid_dst,std::rank<copy_type>::value>::interp(grid_dst,key_dst,x,key_src,interp_pos[Tsrc::value]);
	}
};


#endif /* SRC_GRID_STAGGERED_DIST_GRID_COPY_HPP_ */
