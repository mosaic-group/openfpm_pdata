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

template<typename Decomposition, typename garray>
class Decomposition_encap
{
	Decomposition & dec;
	garray & gd_array;

public:

	Decomposition_encap(Decomposition & dec, garray & gd_array)
	:dec(dec),gd_array(gd_array)
	{}

	Decomposition & internal_dec() const
	{
		return dec;
	}

	/*! \brief Start decomposition
	 *
	 */
	void decompose()
	{
		dec.decompose();

		for(size_t i = 0 ; i < gd_array.size() ; i++)
		{
			Ghost<Decomposition::dims,typename Decomposition::stype> gold = gd_array.get(i).getDecomposition().getGhost();
			gd_array.get(i).getDecomposition() = dec.duplicate(gold);
		}
	}

	/*! \brief Refine the decomposition, available only for ParMetis distribution, for Metis it is a null call
	 *
	 * \param ts number of time step from the previous load balancing
	 *
	 */
	void refine(size_t ts)
	{
		dec.refine();

		for(size_t i = 0 ; i < gd_array.size() ; i++)
		{
			Ghost<Decomposition::dims,typename Decomposition::stype> gold = gd_array.get(i).getDecomposition().getGhost();
			gd_array.get(i).getDecomposition() = dec.duplicate(gold);
		}
	}

	/*! \brief Refine the decomposition, available only for ParMetis distribution, for Metis it is a null call
	 *
	 * \param ts number of time step from the previous load balancing
	 *
	 */
	void redecompose(size_t ts)
	{
		dec.redecompose();

		for(size_t i = 0 ; i < gd_array.size() ; i++)
		{
			Ghost<Decomposition::dims,typename Decomposition::stype> gold = gd_array.get(i).getDecomposition().getGhost();
			gd_array.get(i).getDecomposition() = dec.duplicate(gold);
		}
	}

	auto getDistribution() -> decltype(dec.getDistribution())
	{
		return dec.getDistribution();
	}

	Decomposition_encap<Decomposition,garray> operator=(const Decomposition_encap<Decomposition,garray> & de) const
	{
		for(size_t i = 0 ; i < gd_array.size() ; i++)
		{gd_array.get(i).getDecomposition() = de.gd_array.get(i).getDecomposition();}

		return *this;
	}

	bool write(std::string output) const
	{
		return dec.write(output);
	}
};

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

	//! Boundary conditions of the structure
	periodicity<dim> bc;

	//! array of grids
	//
	openfpm::vector<grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>,
								 HeapMemory,
								 typename memory_traits_lin<grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>>::type,
								 memory_traits_lin,
								 openfpm::grow_policy_identity,STD_VECTOR> gd_array;

	//! Type of structure sub-grid iterator
	typedef decltype(device_grid::type_of_subiterator()) device_sub_it;

	//! Type of structure for the grid iterator
	typedef decltype(device_grid::type_of_iterator()) device_it;

	//! Domain iterator for each distributed grid
	openfpm::vector<grid_dist_iterator<dim,device_grid,device_sub_it,FREE>> git;

	//! Domain and ghost iterator for each distributed grid
	openfpm::vector<grid_dist_iterator<dim,device_grid,device_it,FIXED>> git_g;

	//! Iterator for each distributed grid
	openfpm::vector<grid_dist_iterator_sub<dim,device_grid>> git_sub;

	//! Moving offsets
	openfpm::vector<openfpm::vector<offset_mv<dim>>> mv_off;

	//! background level
	T bck;

	/*! \brief Initialize the others levels
	 *
	 * \param n_grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>lvl number of levels
	 * \param g_sz_lvl grid size on each level
	 *
	 */
	void initialize_other(size_t n_lvl, size_t (& g_sz_lvl)[dim])
	{
		for (size_t i = 0; i < n_lvl - 1 ; i++)
		{
			for (size_t j = 0 ; j < dim ; j++)
			{
				if (bc.bc[j] == NON_PERIODIC)
				{g_sz_lvl[j] = (g_sz_lvl[j]-1)*2 + 1;}
				else
				{g_sz_lvl[j] = g_sz_lvl[j]*2;}
			}

			gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(gd_array.get(0).getDecomposition(),g_sz_lvl,g_int));
			gd_array.last().setBackgroundValue(bck);

			gd_array.last().getDecomposition().free_geo_cell();
			gd_array.last().getDecomposition().getDistribution().destroy_internal_graph();
			gd_array.last().getDecomposition().free_fines();
		}

		recalculate_mvoff();
	}

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
		// set boundary consitions to non periodic

		for (size_t i = 0; i < dim ; i++)
		{bc.bc[i] = NON_PERIODIC;}
	}

	/*! \brief Constructor
	 *
	 * \param domain Simulation domain
	 * \param g ghost extension
	 * \param bc boundary conditions
	 *
	 */
	grid_dist_amr(const Box<dim,St> & domain, const Ghost<dim,long int> & g, periodicity<dim> & bc)
	:domain(domain),g_int(g),bc(bc)
	{
	}

	/*! \brief Initialize the amr grid
	 *
	 * \param dec Decomposition (this parameter is useful in case we want to constrain the AMR to an external decomposition)
	 * \param n_lvl maximum number of levels (0 mean no additional levels)
	 * \param g_sz coarsest grid size on each direction
	 *
	 */
	void initLevels(const Decomposition & dec, size_t n_lvl,const size_t (& g_sz)[dim])
	{
		size_t g_sz_lvl[dim];

		for (size_t i = 0; i < dim ; i++)
		{g_sz_lvl[i] = g_sz[i];}

		// Add the coarse level
		gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(dec,g_sz,g_int));
		gd_array.last().setBackgroundValue(bck);

		initialize_other(n_lvl,g_sz_lvl);
	}

	/*! \brief Initialize the amr grid
	 *
	 * \param dec Decomposition (this parameter is useful in case we want to constrain the AMR to an external decomposition)
	 * \param n_lvl maximum number of levels (0 mean no additional levels)
	 * \param g_sz coarsest grid size on each direction
	 *
	 */
	template<typename TT> void initLevels(const Decomposition_encap<Decomposition,TT> & dec, size_t n_lvl,const size_t (& g_sz)[dim])
	{
		initLevels(dec.internal_dec(),n_lvl,g_sz);
	}

	/*! \brief Recalculate the offset array for the moveLvlUp and moveLvlDw
	 *
	 *
	 *
	 */
	void recalculate_mvoff()
	{
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

	/*! \brief Initialize the amr grid
	 *
	 * \param n_lvl maximum number of levels (0 mean no additional levels)
	 * \param g_sz coarsest grid size on each direction
	 * \param opt options
	 *
	 */
	void initLevels(size_t n_lvl,const size_t (& g_sz)[dim], size_t opt = 0)
	{
		size_t g_sz_lvl[dim];

		for (size_t i = 0; i < dim ; i++)
		{g_sz_lvl[i] = g_sz[i];}

		// Add the coarse level
		gd_array.add(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g_sz,domain,g_int,bc,opt));

		initialize_other(n_lvl,g_sz_lvl);
	}

	/*! \brief Add the computation cost on the decomposition using a resolution function
	 *
	 *
	 * \param md Model to use
	 * \param ts It is an optional parameter approximately should be the number of ghost get between two
	 *           rebalancing at first decomposition this number can be ignored (default = 1) because not used
	 *
	 */
	template <typename Model>inline void addComputationCosts(Model md=Model(), size_t ts = 1)
	{
		gd_array.get(0).addComputationCosts(md,ts);
	}

	/*! \brief Get the object that store the information about the decomposition
	 *
	 * \return the decomposition object
	 *
	 */
	Decomposition_encap<Decomposition,decltype(gd_array)> getDecomposition()
	{
		Decomposition_encap<Decomposition,decltype(gd_array)> tmp(gd_array.get(0).getDecomposition(),gd_array);

		return tmp;
	}

	/*! \brief Get the underlying grid level
	 *
	 * \param lvl level
	 *
	 * \return the grid level
	 *
	 */
	grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> & getLevel(size_t lvl)
	{
		return gd_array.get(lvl);
	}


	grid_dist_amr_key_iterator<dim,device_grid,
							   decltype(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>::type_of_subiterator()),
							   decltype(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>::type_of_subiterator()) >
	getDomainIteratorCells()
	{
		git_sub.clear();

		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			grid_key_dx<dim> start;
			grid_key_dx<dim> stop;

			for (size_t j = 0 ; j < dim ; j++)
			{
				start.set_d(j,0);
				if (bc.bc[j] == NON_PERIODIC)
				{stop.set_d(j,getGridInfoVoid(i).size(j) - 2);}
				else
				{stop.set_d(j,getGridInfoVoid(i).size(j) - 1);}
			}

			git_sub.add(gd_array.get(i).getSubDomainIterator(start,stop));
		}

		return grid_dist_amr_key_iterator<dim,device_grid,
											decltype(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>::type_of_subiterator()),
											decltype(grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>::type_of_subiterator())>(git_sub);
	}

	grid_dist_iterator_sub<dim,device_grid> getDomainIteratorCells(size_t lvl)
	{
		grid_key_dx<dim> start;
		grid_key_dx<dim> stop;

		for (size_t j = 0 ; j < dim ; j++)
		{
			start.set_d(j,0);
			if (bc.bc[j] == NON_PERIODIC)
			{stop.set_d(j,getGridInfoVoid(lvl).size(j) - 2);}
			else
			{stop.set_d(j,getGridInfoVoid(lvl).size(j) - 1);}
		}

		return gd_array.get(lvl).getSubDomainIterator(start,stop);
	}

	/*! \brief Get an iterator to the grid
	 *
	 * \return an iterator to the grid
	 *
	 */
	auto getGridGhostIterator(size_t lvl) -> decltype(gd_array.get(lvl).getGridGhostIterator(grid_key_dx<dim>(),grid_key_dx<dim>()))
	{
		grid_key_dx<dim> key_start;
		grid_key_dx<dim> key_stop;

		for (size_t i = 0 ; i < dim ; i++)
		{
			key_start.set_d(i,g_int.getLow(i));
			key_stop.set_d(i,g_int.getHigh(i) + getGridInfoVoid(lvl).size(i) -1);
		}

		return gd_array.get(lvl).getGridGhostIterator(key_start,key_stop);
	}

	/*! \brief Get an iterator to the grid
	 *
	 * \return an iterator to the grid
	 *
	 */
	auto getGridIterator(size_t lvl) -> decltype(gd_array.get(lvl).getGridIterator())
	{
		return gd_array.get(lvl).getGridIterator();
	}

	/*! \brief Get an iterator to the grid
	 *
	 * \return an iterator to the grid
	 *
	 */
	auto getGridIteratorCells(size_t lvl) -> decltype(gd_array.get(lvl).getGridIterator())
	{
		grid_key_dx<dim> start;
		grid_key_dx<dim> stop;

		for (size_t j = 0 ; j < dim ; j++)
		{
			start.set_d(j,0);
			if (bc.bc[j] == NON_PERIODIC)
			{stop.set_d(j,getGridInfoVoid(lvl).size(j) - 2);}
			else
			{stop.set_d(j,getGridInfoVoid(lvl).size(j) - 1);}
		}

		return gd_array.get(lvl).getGridIterator(start,stop);
	}


	/*! \brief return an iterator over the level lvl
	 *
	 * \param lvl level
	 *
	 * \return an iterator over the level lvl selected
	 *
	 */
	grid_dist_iterator<dim,device_grid,decltype(device_grid::type_of_subiterator()),FREE>
	getDomainIterator(size_t lvl) const
	{
		return gd_array.get(lvl).getDomainIterator();
	}

    /*! \brief return an iterator over the level lvl
     *
     * \param lvl level
     *
     * \return an iterator over the level lvl selected
     *
     */
    grid_dist_iterator<dim,device_grid,
	decltype(device_grid::type_of_iterator()),
    FIXED>
    getDomainGhostIterator(size_t lvl) const
    {
            return gd_array.get(lvl).getDomainGhostIterator();
    }

	/*! \brief Get domain iterator
	 *
	 * \return an iterator over all the grid levels
	 *
	 */
	grid_dist_amr_key_iterator<dim,device_grid, decltype(device_grid::type_of_subiterator())>
	getDomainIterator()
	{
		git.clear();

		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			git.add(gd_array.get(i).getDomainIterator());
		}

		return grid_dist_amr_key_iterator<dim,device_grid,decltype(device_grid::type_of_subiterator())>(git);
	}

    /*! \brief Get domain iterator
     *
     * \return an iterator over all the grid levels
     *
     */
    grid_dist_amr_key_iterator<dim,device_grid, decltype(device_grid::type_of_iterator()),
    		                   grid_dist_iterator<dim,device_grid,decltype(device_grid::type_of_iterator()),FIXED>>
    getDomainGhostIterator()
    {
            git_g.clear();

            for (size_t i = 0 ; i < gd_array.size() ; i++)
            {
                    git_g.add(gd_array.get(i).getDomainGhostIterator());
            }

            return grid_dist_amr_key_iterator<dim,device_grid,decltype(device_grid::type_of_iterator()),
            		                          grid_dist_iterator<dim,device_grid,decltype(device_grid::type_of_iterator()),FIXED>>(git_g);
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


	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto get(size_t lvl, const grid_dist_key_dx<dim> & v1) const -> typename std::add_lvalue_reference<decltype(gd_array.get(lvl).template get<p>(v1))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gd_array.get(lvl).template get<p>(v1);
	}

	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto get(size_t lvl, const grid_dist_key_dx<dim> & v1) -> typename std::add_lvalue_reference<decltype(gd_array.get(lvl).template get<p>(v1))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gd_array.get(lvl).template get<p>(v1);
	}

	//////////////////// Insert functions


	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>
	inline auto insert(const grid_dist_amr_key<dim> & v1)
	-> typename std::add_lvalue_reference<decltype(gd_array.get(v1.getLvl()).template insert<p>(v1.getKey()))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gd_array.get(v1.getLvl()).template insert<p>(v1.getKey());
	}



	/*! \brief Get the reference of the selected element
	 *
	 * \tparam p property to get (is an integer)
	 * \param v1 grid_key that identify the element in the grid
	 *
	 * \return the selected element
	 *
	 */
	template <unsigned int p>inline auto insert(size_t lvl, const grid_dist_key_dx<dim> & v1)
	-> typename std::add_lvalue_reference<decltype(gd_array.get(lvl).template insert<p>(v1))>::type
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		return gd_array.get(lvl).template insert<p>(v1);
	}

	//////////////////////////////////////

	/*! \brief Get the internal distributed grid
	 *
	 * \param lvl level
	 *
	 * \return the internal distributed grid
	 *
	 */
	grid_dist_id<dim,St,T,Decomposition,Memory,device_grid> & getDistGrid(size_t lvl)
	{
		return gd_array.get(lvl);
	}

	//////////////////// Remove functions

	/*! \brief Remove a grid point (this function make sense only in case of
	 *         sparse grid)
	 *
	 * \param v1 grid_key that identify the element in the AMR grid to eleminate
	 *
	 */
	inline void remove(const grid_dist_amr_key<dim> & v1)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		gd_array.get(v1.getLvl()).remove(v1.getKey());
	}

	/*! \brief Remove a grid point (this function make sense only in case of
	 *         sparse grid)
	 *
	 * \param v1 grid_key that identify the element in the AMR grid to eleminate
	 *
	 */
	void remove(size_t lvl, const grid_dist_key_dx<dim> & v1)
	{
#ifdef SE_CLASS2
		check_valid(this,8);
#endif
		gd_array.get(lvl).remove(v1);
	}

	//////////////////////////////////////

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

	/*! \brief It move all the grid parts that do not belong to the local processor to the respective processor
	 *
	 */
	void map(size_t opt = 0)
	{
		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			gd_array.get(i).map();
		}

		recalculate_mvoff();
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

	/*! \brief Return the number of inserted points on a particular level
	 *
	 * \return the number of inserted points
	 *
	 */
	size_t size_inserted(size_t lvl)
	{
		return gd_array.get(lvl).size_local_inserted();
	}

	/*! \brief set the background value
	 *
	 * You can use this function make sense in case of sparse in case of dense
	 * it does nothing
	 *
	 */
	void setBackgroundValue(T & bv)
	{
		for (size_t i = 0 ; i < getNLvl() ; i++)
		{gd_array.get(i).setBackgroundValue(bv);}

		meta_copy<T>::meta_copy_(bv,bck);
	}

	/*! \brief delete all the points in the grid
	 *
	 * In case of sparse grid in delete all the inserted points, in case
	 * of dense it does nothing
	 *
	 */
	void clear()
	{
		for (size_t i = 0 ; i < getNLvl() ; i++)
		{gd_array.get(i).clear();}
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

	/*! \brief Get an object containing the grid informations for a specific level
	 *
	 * \param lvl level
	 *
	 * \return an information object about this grid
	 *
	 */
	const size_t size(size_t lvl,size_t i) const
	{
		return gd_array.get(lvl).size(i);
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

	/*! \brief From a distributed key it return a AMR key that contain also the grid level
	 *
	 * \param lvl level
	 * \param key distributed key
	 *
	 */
	inline grid_dist_amr_key<dim> getAMRKey(size_t lvl, grid_dist_key_dx<dim> key)
	{
		return grid_dist_amr_key<dim>(lvl,key);
	}

	/*! \brief Move up (to coarser level) the key
	 *
	 * \param key multi-resolution AMR key
	 *
	 */
	void moveLvlUp(grid_dist_amr_key<dim> & key)
	{
#ifdef SE_CLASS1

		if (key.getLvl() == 0)
		{std::cerr << __FILE__ << ":" << __LINE__ << " error: we are already at the top level, we cannot go one level up" << std::endl;}

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
	bool write(std::string output, size_t opt = VTK_WRITER | FORMAT_ASCII )
	{
		bool ret = true;

		for (size_t i = 0 ; i < gd_array.size() ; i++)
		{
			ret &= gd_array.get(i).write(output + "_" + std::to_string(i),opt);
		}

		return ret;
	}
};

template<unsigned int dim, typename St, typename T>
using sgrid_dist_amr = grid_dist_amr<dim,St,T,AMR_IMPL_TRIVIAL,CartDecomposition<dim,St>,HeapMemory,sgrid_cpu<dim,T,HeapMemory>>;

#endif /* AMR_GRID_AMR_DIST_HPP_ */
