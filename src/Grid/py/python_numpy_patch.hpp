/*
 * python_patch.hpp
 *
 *  Created on: Sep 5, 2020
 *      Author: i-bird
 */

#ifndef PYTHON_PATCH_HPP_
#define PYTHON_PATCH_HPP_

#include "config.h"

#ifdef HAVE_PYTHON_SERVER

#include <Python.h>
#include <numpy/numpyconfig.h>
#include <numpy/arrayobject.h>
#include <type_traits>

template<typename tpy>
int getPyType()
{
	if (std::is_same<tpy,float>::value)
	{return NPY_FLOAT;}
	else if (std::is_same<tpy,double>::value)
	{return NPY_DOUBLE;}

	return NPY_NOTYPE;
}

template<unsigned int n_prop>
struct py_obj_arr
{
	PyObject * npy[n_prop];
};

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to create a numpy PyObject for every buffer
 *
 *
 */
template<typename dev_grid>
struct create_pyobj
{
	//! wrapped numpy patches
	py_obj_arr<dev_grid::value_type::max_prop> & numpy_patch;

	//! patches
	dev_grid & patch;

	/*! \brief constructor
	 *
	 * \param numpy_patches
	 * \param dst source encapsulated object
	 *
	 */
	inline create_pyobj(py_obj_arr<dev_grid::value_type::max_prop> & numpy_patch,
			            dev_grid & patch)
	:numpy_patch(numpy_patch),patch(patch)
	{
	};

	//! It call the copy function for each property
	template<typename T>
	inline void operator()(T& t) const
	{
		typedef boost::mpl::at<typename dev_grid::value_type::type,T> aggr_sel;

		npy_intp dims[dev_grid::dims+std::rank<typename aggr_sel::type>::value];

		for (int j = 0 ; j < dev_grid::dims ; j++)
		{dims[j] = patch.getGrid().size(j);}

		for (int j = 0 ; j < std::rank<typename aggr_sel::type>::value ; j++)
		{
			if (j == 0)
			{dims[j+dev_grid::dims] = std::extent<typename aggr_sel::type, 0>::value;}

			if (j == 1)
			{dims[j+dev_grid::dims] = std::extent<typename aggr_sel::type, 1>::value;}

			if (j == 2)
			{dims[j+dev_grid::dims] = std::extent<typename aggr_sel::type, 2>::value;}
		}

		int py_type = getPyType<typename std::remove_all_extents<typename aggr_sel::type>::type>();

		PyObject * c = PyArray_SimpleNewFromData(dev_grid::dims,(npy_intp *)&dims,py_type,(void *)patch.template getPointer<T::value>());

		Py_INCREF(c);
	}
};


template<typename dev_grid>
class numpy_patch_class
{
	openfpm::vector<py_obj_arr<dev_grid::value_type::max_prop>> numpy_patches;

public:

	template<typename grid_type>
	inline int create_numpy_patches(openfpm::vector<dev_grid> & dev_g, grid_type & g)
	{
		// This is a macro it return error in case it fail;
		if(PyArray_API == NULL)
		{import_array();}

		// notify
		structs_list sl;
		sl.dim = dev_grid::dims;

		for (int i = 0 ; i < dev_grid::dims ; i++)
		{sl.sizes.add(g.size(i));}

		sl.name = std::string("grid_dist_id");
		sl.type = std::string("native_" + std::string(demangle(typeid(typename grid_type::memory_type).name())) + "_" +
				              std::string(demangle(typeid(typename grid_type::d_grid).name())));

		add_structure(sl);

		numpy_patches.resize(dev_g.size());

		for (int i = 0 ; i < dev_g.size() ; i++)
		{
			create_pyobj<dev_grid> cpy(numpy_patches.get(i),dev_g.get(i));

			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,dev_grid::value_type::max_prop>>(cpy);
		}
	}
};



#else

/*! \brief this class is a functor for "for_each" algorithm
 *
 * This class is a functor for "for_each" algorithm. For each
 * element of the boost::vector the operator() is called.
 * Is mainly used to create a numpy PyObject for every buffer
 *
 *
 */
template<typename dev_grid>
class numpy_patch_class
{
public:
	template<typename grid_type>
	inline void create_numpy_patches(openfpm::vector<dev_grid> & dev_g, grid_type & g)
	{}
};

#endif

#endif /* PYTHON_PATCH_HPP_ */
