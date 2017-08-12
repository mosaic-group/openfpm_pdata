/*
 * staggered_grid.hpp
 *
 *  Created on: Aug 19, 2015
 *      Author: i-bird
 */

#ifndef SRC_GRID_STAGGERED_DIST_GRID_HPP_
#define SRC_GRID_STAGGERED_DIST_GRID_HPP_

#include "Grid/grid_dist_id.hpp"
#include "staggered_dist_grid_util.hpp"
#include "VTKWriter/VTKWriter.hpp"
#include "staggered_dist_grid_copy.hpp"
#include <boost/mpl/vector_c.hpp>

/*! \brief Implementation of the staggered grid
 *
 * \param dim Dimensionality of the staggered grid
 * \param ele elements object on each dimensional objects, must be a stag_elements
 *
 * \verbatim

		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+
		|     |     |     |     |     |     |
		#  *  #  *  #  *  #  *  #  *  #  *  #
		|     |     |     |     |     |     |
		+--#--+--#--+--#--+--#--+--#--+--#--+

\endverbatim

		In the case of a 2D staggered grid we have 3 (in general dim+1 ) elements

		+ = vertex
		# = edge
		* = volume

        ele = stag_ele<scalar<float>,Point_test<float>,scalar<float>>

        It place a scalar on (*) an object Point_test<float> on (#) and an object scalar<float> on (+)

 *
 *
 *
 */
template<unsigned int dim, typename St, typename T, typename Decomposition=CartDecomposition<dim,St>,typename Memory=HeapMemory , typename device_grid=grid_cpu<dim,T>>
class staggered_grid_dist : public grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>
{
	//! position of the properties in the grid cell
	openfpm::vector<comb<dim>> c_prp[T::max_prop_real];

public:

	//! Properties for each grid point
	typedef T value_type;

	//! Number of dimensions
	static const unsigned int dims = dim;

	/*! \brief This constructor is special, it construct an expanded grid that perfectly overlap with the previous
	 *
	 * The key-word here is "perfectly overlap". Using the default constructor you could create
	 * something similar, but because of rounding-off error it can happen that it is not perfectly overlapping
	 *
	 * \param g previous grid
	 * \param gh Ghost part in grid units
	 * \param ext extension of the grid (must be positive on every direction)
	 *
	 */
	template<typename H> staggered_grid_dist(const grid_dist_id<dim,St,H,typename Decomposition::base_type,Memory,grid_cpu<dim,H>> & g,
			                          const Ghost<dim,long int> & gh,
									  Box<dim,size_t> ext)
	:grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g,gh,ext)
	{
	}

	/*! \brief Constructor
	 *
	 * \param g_sz size of the staggered grid
	 * \param domain domain
	 * \param ghost part
	 *
	 *
	 */
	staggered_grid_dist(const size_t (& g_sz)[dim],
			            const Box<dim,St> & domain,
						const Ghost<dim,St> & ghost)
	:grid_dist_id<dim,St,T,Decomposition,Memory,device_grid>(g_sz,domain,ghost)
	{}

	/*! \brief set the staggered positions of the properties
	 *
	 * \tparam property p
	 *
	 * \param cmb a vector containing for each component the position in the cell-grid
	 *
	 */
	template<unsigned int p> void setStagPosition(openfpm::vector<comb<dim>> & cmb)
	{
#ifdef SE_CLASS1
		if (extends< typename boost::mpl::at<typename T::type,boost::mpl::int_<p> >::type >::mul() != cmb.size())
			std::cerr << __FILE__ << ":" << __LINE__ << " error properties has " << extends< typename boost::mpl::at<typename T::type,boost::mpl::int_<p> >::type >::mul() << " components, but " << cmb.size() << "has been defined \n";
#endif
		c_prp.get(p) = cmb;
	}

	/*! \brief It set all the properties defined to be staggered on the default location
	 *
	 */
	void setDefaultStagPosition()
	{
		// for each properties

		stag_set_position<dim,typename T::type_real> ssp(c_prp);

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop_real> >(ssp);
	}

	/*! \brief Copy the staggered grid into a normal one
	 *
	 *
	 * \tparam Grid_dst type of the destination Grid
	 * \tparam pos destination grid properties to fill
	 *
	 * \param g_dst destination grid
	 * \param pd padding of the grid compared to the destination grid
	 * \param start starting point
	 * \param stop stop point
	 *
	 */
	template<typename Grid_dst ,unsigned int ... pos> bool to_normal(Grid_dst & g_dst, const Padding<dim> & pd, const long int (& start)[dim], const long int (& stop)[dim])
	{
		// interpolation points for each property
		openfpm::vector<std::vector<comb<dim>>> interp_pos[sizeof...(pos)];

		typedef boost::mpl::vector_c<unsigned int,pos ... > v_pos_type;

		interp_points<dim,v_pos_type,typename Grid_dst::value_type::type> itp(interp_pos,c_prp);
		boost::mpl::for_each_ref<v_pos_type>(itp);

		// shift the start and stop by the padding
		grid_key_dx<dim> start_k = grid_key_dx<dim>(start);
		grid_key_dx<dim> stop_k = grid_key_dx<dim>(stop);

		// sub-grid iterator over the grid map
		auto g_map_it = this->getSubDomainIterator(start_k,stop_k);

		// Iterator over the destination grid
		auto g_dst_it = g_dst.getDomainIterator();

		// Check that the 2 iterator has the same size
		checkIterator<Grid_dst,decltype(g_map_it),decltype(g_dst_it)>(g_map_it,g_dst_it);

		while (g_map_it.isNext() == true)
		{
			typedef typename to_boost_vmpl<pos...>::type vid;
			typedef boost::mpl::size<vid> v_size;

			auto key_src = g_map_it.get();

			// destination point
			auto key_dst = g_dst_it.get();

			// Transform this id into an id for the Eigen vector

			interp_ele<vid,Grid_dst,typename std::remove_reference<decltype(*this)>::type,sizeof...(pos)> cp(key_dst,g_dst,*this,key_src,interp_pos);

			// For each property in the target grid
			boost::mpl::for_each_ref<boost::mpl::range_c<int,0,v_size::value>>(cp);

			++g_map_it;
			++g_dst_it;
		}

		return true;
	}

	/*! \brief Get the staggered positions
	 *
	 * \return for each property it contain a vector that specify where in the cell
	 *         grid the properties live
	 *
	 */
	const openfpm::vector<comb<dim>>  (& getStagPositions()) [T::max_prop]
	{
		return c_prp;
	}

	/*! \brief Write a vtk file with the information of the staggered grid
	 *
	 * \param str vtk output file
	 *
	 */
	void write(std::string str)
	{
		stag_create_and_add_grid<dim,staggered_grid_dist<dim,St,T,Decomposition,Memory,device_grid>,St> sgw(*this, this->getVC().getProcessUnitID());

		boost::mpl::for_each_ref< boost::mpl::range_c<int,0,T::max_prop> >(sgw);
	}

	/*! \brief Return if the properties is a staggered property or not
	 *
	 * \param prp property to check
	 *
	 * \return true if the property is staggered
	 *
	 */
	bool is_staggered_prop(size_t prp)
	{
		return c_prp[prp].size() != 0;
	}

	/*! \brief Return if the grid is staggered
	 *
	 * \return true
	 *
	 */
	bool is_staggered()
	{
		return true;
	}

	friend class stag_create_and_add_grid<dim,staggered_grid_dist<dim,St,T,Decomposition,Memory,device_grid>,St>;
};

#endif /* SRC_GRID_STAGGERED_DIST_GRID_HPP_ */
