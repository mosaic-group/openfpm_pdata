/*
 * staggered_util.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_DIST_GRID_UTIL_HPP_
#define SRC_GRID_STAGGERED_DIST_GRID_UTIL_HPP_

#include "util/common.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "util/convert.hpp"


/*! \brief write a property that has attributes
 *
 * \tparam ele object we are writing
 * \tparam vtk vtk writer
 * \tparam true in case the basic object has attributes
 *
 */
template<typename ele, typename vtk, bool has_attributes=has_attributes<ele>::value>
struct vtk_write
{
	/*! \brief Add the grid with attributes name
	 *
	 * \param vtk_w VTK writer
	 * \param output where to write
	 * \param i property to write
	 *
	 */
	vtk_write(vtk vtk_w, const std::string output, const size_t i)
	{
		vtk_w.write(output + "_" + ele::attributes::name[i] + ".vtk",ele::attributes::name[i]);
	}
};

/*! \brief Add to the vtk writer the key
 *
 * \tparam ele object we are writing
 * \tparam vtk vtk writer
 * \tparam false in case the basic object has not attributes
 *
 */
template<typename ele, typename vtk>
struct vtk_write<ele,vtk,false>
{
	/*! \brief Add the grid with attributes name
	 *
	 * \param vtk_w VTK writer
	 * \param output where to write
	 * \param i property to write
	 *
	 */
	vtk_write(vtk vtk_w, const std::string output, const size_t i)
	{
		vtk_w.write(output + "_" + std::to_string(i) + ".vtk","attr" + std::to_string(i));
	}
};


/*! \brief Classes to get the number of components of the properties
 *
 */
template<typename T>
struct extends
{
	/*! \brief Scalar case
	 *
	 * \return 1 component
	 *
	 */
	static inline size_t mul()
	{
		return 1;
	}

	/*! \brief Dimensionality
	 *
	 * \return 0
	 *
	 */
	static inline size_t dim()
	{
		return 0;
	}
};

//! Partial specialization for N=1 1D-Array
template<typename T,size_t N1>
struct extends<T[N1]>
{
	/*! \brief Vector case return N1 component
	 *
	 * \return N1
	 *
	 */
	static inline size_t mul()
	{
		return N1;
	}

	/*! Dimensionality 1
	 *
	 * \return 1
	 *
	 */
	static inline size_t dim()
	{
		return 1;
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct extends<T[N1][N2]>
{
	/*! \brief Matrix case return N1*N2 component
	 *
	 * \return N1*N2
	 *
	 */
	static inline size_t mul()
	{
		return N1 * N2;
	}

	/*! Dimensionality 2
	 *
	 * \return 2
	 *
	 */
	static inline size_t dim()
	{
		return 2;
	}
};

//! Partial specialization for N=3
template<typename T,size_t N1,size_t N2,size_t N3>
struct extends<T[N1][N2][N3]>
{
	//! number of elements
	static inline size_t mul()
	{
		return N1 * N2 * N3;
	}

	/*! number of indexes
	 *
	 * \return 3
	 *
	 */
	static inline size_t dim()
	{
		return 3;
	}
};

//! Partial specialization for N=4
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4>
struct extends<T[N1][N2][N3][N4]>
{
	//! number of elements
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4;
	}

	/*! number of indexes
	 *
	 * \return 4
	 *
	 */
	static inline size_t dim()
	{
		return 4;
	}
};

//! Partial specialization for N=5
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5>
struct extends<T[N1][N2][N3][N4][N5]>
{
	//! number of elements
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5;
	}

	/*! number of indexes
	 *
	 * \return 5
	 *
	 */
	static inline size_t dim()
	{
		return 5;
	}
};

//! Partial specialization for N=6
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6>
struct extends<T[N1][N2][N3][N4][N5][N6]>
{
	//! number of elements
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6;
	}

	/*! number of indexes
	 *
	 * \return 6
	 *
	 */
	static inline size_t dim()
	{
		return 6;
	}
};

//! Partial specialization for N=7
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7>
struct extends<T[N1][N2][N3][N4][N5][N6][N7]>
{
	//! number of elements
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7;
	}

	/*! number of indexes
	 *
	 * \return 7
	 *
	 */
	static inline size_t dim()
	{
		return 7;
	}
};

//! Partial specialization for N=8
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8>
struct extends<T[N1][N2][N3][N4][N5][N6][N7][N8]>
{
	/*! number of elements
	 *
	 * \return the number of elements as N1*N2*N3*.........
	 *
	 */
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8;
	}

	/*! number of indexes
	 *
	 * \return 8
	 *
	 */
	static inline size_t dim()
	{
		return 8;
	}
};

//! Partial specialization for N=9
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9>
struct extends<T[N1][N2][N3][N4][N5][N6][N7][N8][N9]>
{
	/*! number of elements
	 *
	 * \return the number of elements as N1*N2*N3*.........
	 *
	 */
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8 * N9;
	}

	/*! number of indexes
	 *
	 * \return 9
	 *
	 */
	static inline size_t dim()
	{
		return 9;
	}
};

//! Partial specialization for N=10
template<typename T,size_t N1,size_t N2,size_t N3,size_t N4,size_t N5, size_t N6, size_t N7, size_t N8, size_t N9, size_t N10>
struct extends<T[N1][N2][N3][N4][N5][N6][N7][N8][N9][N10]>
{
	/*! number of elements
	 *
	 * \return the number of elements as N1*N2*N3*.........
	 *
	 */
	static inline size_t mul()
	{
		return N1 * N2 * N3 * N4 * N5 * N6 * N7 * N8 * N9 * N10;
	}

	/*! number of indexes
	 *
	 * \return 10
	 *
	 */
	static inline size_t dim()
	{
		return 10;
	}
};

///////////////////// Copy grid extends

/*! \brief Classes to copy each component into a grid and add to the VTKWriter the grid
 *
 * \param T property to write
 * \param dim dimansionality
 * \param St type of space
 * \param VTK VTK writer
 *
 */
template<typename T>
struct write_stag
{
	/*! \brief write the staggered grid
	 *
	 * \tparam p_val property we are going to write
	 * \tparam sg staggered grid type
	 * \tparam v_g vector of grids
	 *
	 * \param st_g staggered grid
	 * \param v_g vector of grids
	 * \param lg local grid of the staggered grid we are writing
	 *
	 */
	template<unsigned int p_val, typename sg, typename v_g> static inline void write(sg & st_g, v_g & vg,size_t lg)
	{
		// Add a grid;
		vg.add();
		size_t k = vg.size() - 1;

		// Get the source and destination grid
		auto & g_src = st_g.get_loc_grid(lg);
		auto & g_dst = vg.get(k);

		// Set dimensions and memory
		g_dst.resize(g_src.getGrid().getSize());

		// copy

		auto it = vg.get(k).getIterator();

		while(it.isNext())
		{
			g_dst.template get<0>(it.get()) = g_src.template get<p_val>(it.get());

			++it;
		}
	}
};

/*! \brief for each component add a grid fill it, and add to the VTK writer
 *
 * \param T Property to copy
 * \param N1 number of components
 *
 */
template<typename T,size_t N1>
struct write_stag<T[N1]>
{
	/*! \brief write the staggered grid
	 *
	 * \tparam p_val property we are going to write
	 * \tparam sg staggered grid type
	 * \tparam v_g vector of grids
	 *
	 * \param st_g staggered grid
	 * \param v_g vector of grids
	 * \param lg local grid of the staggered grid we are writing
	 *
	 */
	template<unsigned int p_val, typename sg, typename v_g> static inline void write(sg & st_g, v_g & vg,size_t lg)
	{
		for (size_t i = 0 ; i < N1 ; i++)
		{
			// Add a grid;
			vg.add();
			size_t k = vg.size() - 1;

			// Get the source and destination grid
			auto & g_src = st_g.get_loc_grid(lg);
			auto & g_dst = vg.get(k);

			// Set dimensions and memory
			g_dst.resize(g_src.getGrid().getSize());

			auto it = vg.get(k).getIterator();

			while(it.isNext())
			{
				g_dst.template get<0>(it.get()) = g_src.template get<p_val>(it.get())[i];

				++it;
			}
		}
	}
};

//! Partial specialization for N=2 2D-Array
template<typename T,size_t N1,size_t N2>
struct write_stag<T[N1][N2]>
{
	/*! \brief write the staggered grid
	 *
	 * \tparam p_val property we are going to write
	 * \tparam sg staggered grid type
	 * \tparam v_g vector of grids
	 *
	 * \param st_g staggered grid
	 * \param v_g vector of grids
	 * \param lg local grid of the staggered grid we are writing
	 *
	 */
	template<unsigned int p_val, typename sg, typename v_g> static inline void write(sg & st_g, v_g & vg,size_t lg)
	{
		for (size_t i = 0 ; i < N1 ; i++)
		{
			for (size_t j = 0 ; j < N2 ; j++)
			{
				// Add a grid;
				vg.add();
				size_t k = vg.size() - 1;

				// Set dimensions and memory
				vg.get(k).resize(st_g.get_loc_grid(lg).getGrid().getSize());

				// copy
				auto & g_src = st_g.get_loc_grid(lg);
				auto & g_dst = vg.get(k);
				auto it = vg.get(k).getIterator();

				while(it.isNext())
				{
					g_dst.template get<0>(it.get()) = g_src.template get<p_val>(it.get())[i][j];

					++it;
				}
			}
		}
	}
};

///////////////////// Staggered default positioning ////////////////////////

/*! \brief this class is a functor for "for_each" algorithm
 *
 * For each element of the boost::vector the operator() is called.
 * Is mainly used to produce a default position vector for each
 * property
 *
 * \tparam dim dimensionality
 * \tparam v boost::fusion::vector of properties
 * \tparam has_posMask case when v has a position mask
 *
 */

template<unsigned int dim, typename v, bool has_pM = has_posMask<v>::value>
class stag_set_position
{
	//! vector containing the position of the properties in the cells (staggered properties are staggered)
	// within the cell
	openfpm::vector<comb<dim>> (& pos_prp)[boost::fusion::result_of::size<v>::type::value];

public:

	/*! \brief Constructor
	 *
	 * \param vector of the staggered position (It is going to be filled by this class)
	 *
	 */
	stag_set_position( openfpm::vector<comb<dim>> (& pos_prp)[boost::fusion::result_of::size<v>::type::value])
	:pos_prp(pos_prp)
	{}

	/*! It calculate the staggered position for every property
	 *
	 * \param t property
	 *
	 */
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
				c.value(i) = -1;

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
				comb<dim> c1 = pos_prp[T::value-1].get(i);
				for (size_t j = 0 ; j < dim ; j++)
				{
					comb<dim> c2;
					c2.zero();
					c2.value(i) = 1;

					comb<dim> c_res = (c1 + c2) & 0x1;

					pos_prp[T::value].add(c_res);
				}
			}
		}
		else if (dim_prp > 2)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Tensor of rank bigger than 2 are not supported";
		}
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

template<unsigned int dim, typename v>
class stag_set_position<dim,v,false>
{
private:

	//! vector containing the position of the properties in the cells (staggered properties are staggered)
	// within the cell
	openfpm::vector<comb<dim>> (& pos_prp)[boost::fusion::result_of::size<v>::type::value];


public:

	/*! \brief Constructor
	 *
	 * \param vector of the staggered position (It is going to be filled by this class)
	 *
	 */
	stag_set_position( openfpm::vector<comb<dim>> (& pos_prp)[boost::fusion::result_of::size<v>::type::value])
	:pos_prp(pos_prp)
	{}

	/*! It calculate the staggered position for every property
	 *
	 * \param t property
	 *
	 */
	template<typename T>
	void operator()(T& t) const
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at<v,typename boost::mpl::int_<T::value>>::type prop;

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
				c.getComb()[i] = -1;

				pos_prp[T::value].add(c);
			}
		}
		else if (dim_prp == 2)
		{
			// Diagonal part live in
			for (size_t i = 0 ; i < dim ; i++)
			{
				comb<dim> c1 = pos_prp[T::value-1].get(i);
				for (size_t j = 0 ; j < dim ; j++)
				{
					comb<dim> c2;
					c2.zero();
					c2.getComb()[j] = 1;

					comb<dim> c_res = (c2 + c1).flip();

					pos_prp[T::value].add(c_res);
				}
			}
		}
		else if (dim_prp > 2)
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " Tensor of rank bigger than 2 are not supported";
		}
	}
};

/*! \brief It create separated grid for each properties to write them into a file
 *
 * \tparam dim dimensionality of the grids
 * \tparam obj type object to print, must be in OpenFPM format
 *
 */
template<unsigned int dim, typename st_grid, typename St>
class stag_create_and_add_grid
{

	size_t p_id;

	// staggered grid to write
	st_grid & st_g;

public:

	/*! \brief Constructor
	 *
	 * \param st_g staggered grid
	 * \param p_id process id
	 *
	 */
	stag_create_and_add_grid(st_grid & st_g, size_t p_id)
	:p_id(p_id),st_g(st_g)
	{}

	template<unsigned int p_val> void out_normal()
	{
		// property type
		typedef typename boost::mpl::at< typename st_grid::value_type::type , typename boost::mpl::int_<p_val> >::type ele;

		// create an openfpm format object from the property type
		typedef object<typename boost::fusion::vector<ele>> d_object;

		VTKWriter<boost::mpl::pair<grid_cpu<dim, d_object >,St>,VECTOR_GRIDS> vtk_w;

		// Create a vector of grids

		openfpm::vector< grid_cpu<dim, d_object > > vg(st_g.getN_loc_grid());

		// for each domain grid
		for (size_t i = 0 ; i < vg.size() ; i++)
		{
			// Set dimansions and memory
			vg.get(i).resize(st_g.get_loc_grid(i).getGrid().getSize());

			auto & g_src = st_g.get_loc_grid(i);
			auto & g_dst = vg.get(i);

			auto it = vg.get(i).getIterator();

			while(it.isNext())
			{
				object_si_d< decltype(g_src.get_o(it.get())),decltype(g_dst.get_o(it.get())) ,OBJ_ENCAP,p_val>(g_src.get_o(it.get()),g_dst.get_o(it.get()));

				++it;
			}

			Point<dim,St> offset = st_g.getOffset(i);
			Point<dim,St> spacing = st_g.getSpacing();
			Box<dim,size_t> dom = st_g.getDomain(i);

			vtk_w.add(g_dst,offset,spacing,dom);
		}

		vtk_w.write("vtk_grids_st_" + std::to_string(p_id) + "_" + std::to_string(p_val) + ".vtk",st_g.getPropNames(),"grids",file_type::BINARY);
	}

	template<unsigned int p_val> void out_staggered()
	{
		// property type
		typedef typename boost::mpl::at< typename st_grid::value_type::type , typename boost::mpl::int_<p_val> >::type ele;

		// Eliminate the extends
		typedef typename std::remove_all_extents<ele>::type r_ele;

		// create an openfpm format object from the property type
		typedef object<typename boost::fusion::vector<r_ele>> d_object;

		VTKWriter<boost::mpl::pair<grid_cpu<dim, d_object >,St>,VECTOR_ST_GRIDS> vtk_w;

		// Create a vector of grids
		openfpm::vector< grid_cpu<dim, d_object > > vg;
		vg.reserve(st_g.getN_loc_grid() * extends<ele>::mul());

		size_t k = 0;

		// for each domain grid
		for (size_t i = 0 ; i < st_g.getN_loc_grid() ; i++)
		{
			write_stag<ele>::template write<p_val, st_grid,openfpm::vector< grid_cpu<dim, d_object > > >(st_g,vg,i);

			// for each component
			for ( ; k < vg.size() ; k++)
			{
				Point<dim,St> offset = st_g.getOffset(i);
				Point<dim,St> spacing = st_g.getSpacing();
				Box<dim,size_t> dom = st_g.getDomain(i);

				vtk_w.add(i,vg.get(k),offset,spacing,dom,st_g.c_prp[p_val].get(k));
			}

			k = vg.size();
		}

		vtk_write<typename st_grid::value_type,VTKWriter<boost::mpl::pair<grid_cpu<dim, d_object >,St>,VECTOR_ST_GRIDS>> v(vtk_w,"vtk_grids_st_" + std::to_string(p_id),p_val);
	}

	//! It call the copy function for each property
	template<typename T>
	void operator()(T& t)
	{
		if (st_g.is_staggered_prop(T::value) == false)
			out_normal<T::value>();
		else
			out_staggered<T::value>();
	}
};

/*! \brief Check that the size of the iterators match
 *
 * It check the the boxes that the sub iterator defines has same dimensions, for example
 * if the first sub-iterator, iterate from (1,1) to (5,3) and the second from (2,2) to (6,4)
 * they match (2,2) to (4,6) they do not match
 *
 * \tparam Grid_map type of the map grid
 * \tparam Grid_dst type of the destination grid
 *
 * \param it1 Iterator1
 * \param it2 Iterator2
 *
 * \return true if they match
 *
 */
template<typename Eqs_sys, typename it1_type, typename it2_type> bool checkIterator(const it1_type & it1, const it2_type & it2)
{
#ifdef SE_CLASS1

	grid_key_dx<Eqs_sys::dims> it1_k = it1.getStop() - it1.getStart();
	grid_key_dx<Eqs_sys::dims> it2_k = it2.getStop() - it2.getStart();

	for (size_t i = 0 ; i < Eqs_sys::dims ; i++)
	{
		if (it1_k.get(i) !=  it2_k.get(i))
		{
			std::cerr << __FILE__ << ":" << __LINE__ << " error src iterator and destination iterator does not match in size\n";
			return false;
		}
	}

	return true;
#else

	return true;

#endif
}

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to calculate the interpolation points for each
 * property in a staggered grid
 *
 * \tparam dim Dimensionality
 * \tparam v_prp_id vector of properties id
 * \tparam v_prp_type vector with the properties
 *
 */
template<unsigned int dim, unsigned int n_prop, typename v_prp_id, typename v_prp_type>
struct interp_points
{
/*#ifdef SE_CLASS3

	// number of properties we are processing
	typedef boost::mpl::size<v_prp_id> v_size_old;

	typedef boost::mpl::int_<v_size_old::value-3> v_size;

#else*/

	// number of properties we are processing
	typedef boost::mpl::size<v_prp_id> v_size;

//#endif

	// interpolation points for each property
	openfpm::vector<std::vector<comb<dim>>> (& interp_pts)[v_size::value];

	// staggered position for each property
	const openfpm::vector<comb<dim>> (&stag_pos)[n_prop];

	/*! \brief constructor
	 *
	 * It define the copy parameters.
	 *
	 * \param inter_pts array that for each property contain the interpolation points for each components
	 * \param staggered position for each property and components
	 *
	 */
	inline interp_points(openfpm::vector<std::vector<comb<dim>>> (& interp_pts)[v_size::value],const openfpm::vector<comb<dim>> (&stag_pos)[n_prop])
	:interp_pts(interp_pts),stag_pos(stag_pos){};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t)
	{
		// This is the type of the object we have to copy
		typedef typename boost::mpl::at_c<v_prp_type,T::value>::type prp_type;
		typedef typename boost::mpl::at<v_prp_id,T>::type p_id;

		interp_pts[T::value].resize(stag_pos[p_id::value].size());

		for (size_t i = 0 ; i < stag_pos[p_id::value].size() ; i++)
		{
			// Create the interpolation points
			interp_pts[T::value].get(i) = SubHyperCube<dim,dim - std::rank<prp_type>::value>::getCombinations_R(stag_pos[p_id::value].get(i),0);

			// interp_point are -1,0,1, map the -1 to 0 and 1 to -1
			for (size_t j = 0 ; j < interp_pts[T::value].get(i).size() ; j++)
			{
				for (size_t k = 0 ; k < dim ; k++)
					interp_pts[T::value].get(i)[j].getComb()[k] = - ((interp_pts[T::value].get(i)[j].getComb()[k] == -1)?0:interp_pts[T::value].get(i)[j].getComb()[k]);
			}
		}
	}
};

#endif /* SRC_GRID_STAGGERED_DIST_GRID_UTIL_HPP_ */
