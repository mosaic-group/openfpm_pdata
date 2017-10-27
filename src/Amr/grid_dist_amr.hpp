/*
 * grid_amr_dist.hpp
 *
 *  Created on: Sep 21, 2017
 *      Author: i-bird
 */

#ifndef AMR_GRID_AMR_DIST_HPP_
#define AMR_GRID_AMR_DIST_HPP_

#include "Grid/grid_dist_id.hpp"
#include "Amr/grid_dist_amr_key_iterator.hpp"

#define AMR_IMPL_TRIVIAL 1
#define AMR_IMPL_PATCHES 2
#define AMR_IMPL_OPENVDB 3

template<unsigned int dim, typename St, typename T, unsigned int impl=AMR_IMPL_TRIVIAL ,typename Decomposition = CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T> >
class grid_dist_amr
{

};

/*! \brief It contain the offset necessary to move to coarser and finer level grids
 *
 */
template<unsigned int dim>
struct offset_mv
{
	//! offset to move up on an upper grid (coarse)
	Point<dim,long int> up;

	//! offset to move on the lower grid (finer)
	Point<dim,long int> dw;
};

/*! \brief AMR Adaptive Multi Resolution Grid
 *
 * \tparam dim Dimensionality
 * \tparam St type of space
 * \tparam T what each point of the grid store
 * \tparam Decomposition type of decomposition
 *
 */
template<unsigned int dim, typename St, typename T,typename Decomposition,typename Memory, typename device_grid >
class grid_dist_amr<dim,St,T,AMR_IMPL_TRIVIAL,Decomposition,Memory,device_grid>
{
	//! Simulation domain
	Box<dim,St> domain;

	//! Ghost integer
	Ghost<dim,long int> g_int;

	//! array of grids
	openfpm::vector<grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>> gd_array;

	//! Iterator for each distributed grid
	openfpm::vector<grid_dist_iterator<dim,device_grid,FREE>> git;

	//! Iterator for each distributed grid
	openfpm::vector<grid_dist_iterator_sub<dim,device_grid>> git_sub;

	//! Moving offsets
	openfpm::vector<openfpm::vector<offset_mv<dim>>> mv_off;

public:


	/*! \brief Constructor
	 *
	 * \param domain Simulation domain
	 * \param g ghost extension
	 *
	 */
	grid_dist_amr(const Box<dim,St> & domain, const Ghost<dim,long int> & g)
	:domain(domain),g_int(g)
	{
	}

	/*! \brief Initialize the amr grid
	 *
	 * \param n_lvl maximum number of levels (0 mean no additional levels)
	 * \param g_sz coarsest grid size on each direction
	 *
	 */
	void initLevels(size_t n_lvl,const size_t (& g_sz)[dim])
	{
		size_t g_sz_lvl[dim];

		for (size_t i = 0; i < dim ; i++)
		{g_sz_lvl[i] = g_sz[i];}

		// Add the coarse level
		gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g_sz,domain,g_int));

		for (size_t i = 0; i < n_lvl - 1 ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				g_sz_lvl[j] = (g_sz_lvl[j]-1)*2 + 1;
			}

			gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(gd_array.get(0).getDecomposition(),g_sz_lvl,g_int));
		}

		// Here we calculate the offset to move one level up and one level down
		// in global coordinated moving one level up is multiply the coordinates by 2
		// and moving one level down is dividing by 2. In local coordinates is the same
		// with the exception that because of the decomposition you can have an offset
		// look at the picture below
		//
		// (-1)  (0)
		//  * |   *     *    coarse level
		//  * |*  *  *  *    finer level
		//    |(0)(1)
		//
		//  Line of the decomposition
		//
		// The coarse level point 0 in local coordinates converted to the finer level is not
		// just 2*0 = 0 but is 2*(0) + 1 so a formula like 2*x+offset is required. here we calculate
		// these offset. In the case of moving from finer to coarse is the same the formula is
		// Integer_round(x+1)/2 - 1
		//
		mv_off.resize(gd_array.size());

		for (size_t i = 1 ; i < gd_array.size() ; i++)
		{
			auto & g_box_c = gd_array.get(i-1).getLocalGridsInfo();
			auto & g_box_f = gd_array.get(i).getLocalGridsInfo();

#ifdef SE_CLASS1

			if (g_box_c.size() != g_box_f.size())
			{
				std::cerr << __FILE__ << ":" << __LINE__ << " error it seem that the AMR construction between level " <<
						i << " and " << i-1 << " is inconsistent" << std::endl;
			}

#endif

			mv_off.get(i-1).resize(g_box_f.size());
			mv_off.get(i).resize(g_box_f.size());

			for (size_t j = 0 ; j < g_box_f.size() ; j++)
			{
				for (size_t s = 0 ; s < dim ; s++)
				{
					size_t d_orig_c = g_box_c.get(j).origin.get(s);
					size_t d_orig_f = g_box_f.get(j).origin.get(s);

					mv_off.get(i-1).get(j).dw.get(s) = d_orig_c*2 - d_orig_f;
					mv_off.get(i).get(j).up.get(s) = d_orig_c*2 - d_orig_f;
				}
			}
		}
	}

	/*! \brief Get domain iterator
	 *
	 * \return an iterator over all the grid levels
	 *
	 */
	grid_dist_amr_key_iterator<dim,device_grid> getDomainIterator()
	{
		git.clear();

		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			git.add(gd_array.get(i).getDomainIterator());
		}

		return grid_dist_amr_key_iterator<dim,device_grid>(git);
	}

	/*! \brief return an iterator over the level lvl
	 *
	 * \param lvl level
	 *
	 * \return an iterator over the level lvl selected
	 *
	 */
	grid_dist_iterator<dim,device_grid,FREE> getDomainIterator(size_t lvl) const
	{
		return gd_array.get(lvl).getDomainIterator();
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto get(const grid_dist_amr_key<dim> & v1) const -> typename std::add_lvalue_reference<decltype(gd_array.get(v1.getLvl()).template get<p>(v1.getKey()))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gd_array.get(v1.getLvl()).template get<p>(v1.getKey());
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto get(const grid_dist_amr_key<dim> & v1) -> typename std::add_lvalue_reference<decltype(gd_array.get(v1.getLvl()).template get<p>(v1.getKey()))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gd_array.get(v1.getLvl()).template get<p>(v1.getKey());
	}

	/*! \brief It synchronize the ghost parts
	 *
	 * \tparam prp... Properties to synchronize
	 *
	 */
	template<int... prp> void ghost_get()
	{
		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			gd_array.get(i).template ghost_get<prp...>();
		}
	}

	/*! \brief Apply the ghost put
	 *
	 * \tparam prp... Properties to apply ghost put
	 *
	 */
	template<template<typename,typename> class op,int... prp> void ghost_put()
	{
		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			gd_array.get(i).template ghost_put<op,prp...>();
		}
	}

	/*! \brief Get an object containing the grid informations for a specific level
	 *
	 * \param lvl level
	 *
	 * \return an information object about this grid
	 *
	 */
	const grid_sm<dim,void> & getGridInfoVoid(size_t lvl) const
	{
		return gd_array.get(lvl).getGridInfoVoid();
	}

	/*! \brief Return the maximum number of levels in the AMR struct
	 *
	 * \return the number of levels
	 *
	 */
	size_t getNLvl()
	{
		return gd_array.size();
	}

	/*! \brief Move down (to finer level) the key
	 *
	 * \param key multi-resolution AMR key
	 *
	 */
	void moveLvlDw(grid_dist_amr_key<dim> & key)
	{
#ifdef SE_CLASS1

		if (key.getLvl() >= getNLvl() - 1)
		{std::cerr << __FILE__ << ":" << __LINE__ << " error: we are already at the last level, we cannot go one level down" << std::endl;}

#endif

		auto & key_ref = key.getKeyRef().getKeyRef();
		size_t lvl = key.getLvl();

		for (size_t i = 0 ; i < dim ; i++)
		{
			key_ref.set_d(i,(key_ref.get(i) << 1) + mv_off.get(key.getLvl()).get(key.getKeyRef().getSub()).dw.get(i) );
		}

		key.setLvl(lvl+1);
	}

	/*! \brief Move up (to coarser level) the key
	 *
	 * \param key multi-resolution AMR key
	 *
	 */
	void moveLvlUp(grid_dist_amr_key<dim> & key)
	{
#ifdef SE_CLASS1

		if (key.getLvl() >= getNLvl() - 1)
		{std::cerr << __FILE__ << ":" << __LINE__ << " error: we are already at the last level, we cannot go one level down" << std::endl;}

#endif

		auto & key_ref = key.getKeyRef().getKeyRef();
		size_t lvl = key.getLvl();

		for (size_t i = 0 ; i < dim ; i++)
		{
			key_ref.set_d(i,(key_ref.get(i) - mv_off.get(key.getLvl()).get(key.getKeyRef().getSub()).up.get(i)) >> 1);
		}

		key.setLvl(lvl-1);
	}

	/*! \brief Get the position on the grid in global coordinates
	 *
	 * \param v1 amr key
	 *
	 * \return the position in global coordinates
	 *
	 */
	grid_key_dx<dim> getGKey(const grid_dist_amr_key<dim> & v1)
	{
		return gd_array.get(v1.getLvl()).getGKey(v1.getKey());
	}

	/*! \brief return the spacing for the grid in the level lvl
	 *
	 * \param lvl level
	 *
	 * \return return the spacing
	 *
	 */
	Point<dim,St> getSpacing(size_t lvl)
	{
		return gd_array.get(lvl).getSpacing();
	}

	/*! \brief Write on vtk file
	 *
	 * \param output filename output
	 *
	 */
	bool write(std::string output)
	{
		bool ret = true;

		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			ret &= gd_array.get(i).write(output + "_" + std::to_string(i));
		}

		return ret;
	}
};


#endif /* AMR_GRID_AMR_DIST_HPP_ */
