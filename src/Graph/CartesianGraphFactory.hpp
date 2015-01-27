/*
 * CartesianGraphFactory.hpp
 *
 *  Created on: Nov 28, 2014
 *      Author: i-bird
 */

#ifndef CARTESIANGRAPHFACTORY_HPP_
#define CARTESIANGRAPHFACTORY_HPP_

#include "map_vector.hpp"
#include "map_graph.hpp"
#include "grid.hpp"
#include "Space/Shape/Box.hpp"

/*! \brief This class work as a functor
 *
 * For each number in the boost::mpl::vector (for example 3 6) set the properties of the vertex at the
 * specified id (3 6) with pos[d] * spacing[d] with d running from 0 to 1, pos[d] the position id of the vertex
 * spacing the grid spacing
 *
 * Example
 *
 * if we give a grid_key of dimension 2 4x4 the expression "pos[d] * spacing[d]"
 * will assume the value
 *
 * (0.0 0.0) (0.25 0.0) ...... (1.0 0.0)
 * (0.0 0.25)................. (1.0 0.25)
 * ....................................
 * (0.0 1.0).................. (1.0 1.0)
 *
 * and the properties 3 6 will be filled with the numbers 0.0 0.0    .......  1.0 1.0
 * progressively
 *
 * \tparam dim Dimensionality of the cartesian grid
 * \tparam dT type of the domain
 * \tparam G_v vertex type object
 * \tparam v boost::mpl::vector containing all the index to fill
 * \tparam is_stub when is true, produce a trivial operator(),
 *         to use when v is an empty vector to avoid compilation error
 *
 */

template<unsigned int dim, typename dT, typename G_v, typename v, bool is_stub>
class fill_prop
{
	//! Reference to an array containing the spacing
	const dT (& szd)[dim];

	//! grid_key_dx Reference containing the actual position
	grid_key_dx<dim> & gk;

	//! Vertex object to fill
	G_v & g_v;

public:

	//! Fill the object from where to take the properties
	fill_prop(G_v & g_v , const dT (& szd)[dim], grid_key_dx<dim> & gk)
	:szd(szd),gk(gk),g_v(g_v)
	{}

	//! It call the function for each property we want to copy
    template<typename T>
    void operator()(T& t) const
    {
    	typedef typename boost::fusion::result_of::at<v,boost::mpl::int_<T::value>>::type t_val;

    	g_v.template get<t_val::value>() = gk.get(T::value) * szd[T::value];
    }
};

/*! \brief Graph constructor function specialization
 *
 * On C++ partial function specialization is not allowed, so we need a class to do it
 *
 * \see CartesianGraphFactory method construct
 *
 */

template<unsigned int dim, typename Graph, unsigned int se,typename T, unsigned int dim_c, int... pos>
class Graph_constructor_impl
{
public:
	//! Construct cartesian graph
	static Graph construct(std::vector<size_t> sz, Box<dim,T> dom)
	{
#ifdef DEBUG
		//! The size is wrong signal it

		if (sz.size() != dim)
		{std::cerr << "Error this factory has been specialized for catesian grid of dimension " << dim << "\n";}

#endif

		// Calculate the size of the hyper-cubes on each dimension

		T szd[dim];

		for (int i = 0 ; i < dim ; i++)
		{szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid<dim,void> g(sz);

		// Create a graph with the number of vertices equal to the number of
		// grid point

		//! Graph to construct

		Graph gp(g.size());

		/******************
		 *
		 * Create the edges and fill spatial
		 * information properties
		 *
		 ******************/

		//! Construct a key iterator

		grid_key_dx_iterator<dim> k_it(g);

		//! Iterate through all the elements

		while (k_it.isNext())
		{
			grid_key_dx<dim> key = k_it.get();

			// Vertex object
			auto obj = gp.vertex(g.LinId(key));

			// vertex spatial properties functor
			fill_prop<dim,T,decltype(gp.vertex(g.LinId(key))), typename to_boost_mpl<pos...>::type, sizeof...(pos) == 0 > flp(obj,szd,key);

			// fill properties

			boost::mpl::for_each< boost::mpl::range_c<int,0,sizeof...(pos)> >(flp);

			// Get the combinations of dimension d

			for (int d = dim-1 ; d >= dim_c ; d--)
			{
				// create the edges for that dimension

				std::vector<comb<dim>> c = hc.getCombinations_R(d);

				// for each combination calculate a safe linearization and create an edge

				for (int j = 0 ; j < c.size() ; j++)
				{
					// Calculate the element size

					T ele_sz = 0;

					// for each dimension multiply and reduce

					for (int s = 0 ; s < dim ; s++)
					{
						ele_sz += szd[s] * abs(c[j][s]);
					}

					// Calculate the end point vertex id
					// Calculate the start point id

					size_t start_v = g.LinId(key);
					size_t end_v = g.template LinId<CheckExistence>(key,c[j].getComb());

					// Add an edge and set the the edge property to the size of the face (communication weight)
					gp.template addEdge<CheckExistence>(start_v,end_v).template get<se>() = ele_sz;
				}
			}

			// Fill vertex properties



			++k_it;
		}

		return gp;
	}
};

/*! \brief Graph constructor function specialization
 *
 * On C++ partial function specialization is not allowed, so we need a class to do it
 * This specialization handle the case when we have NO_EDGE option active
 *
 * \see CartesianGraphFactory method construct
 *
 */

template<unsigned int dim, typename Graph,typename T, unsigned int dim_c, int... pos>
class Graph_constructor_impl<dim,Graph,NO_EDGE,T,dim_c,pos...>
{
public:
	//! Construct cartesian graph
	static Graph construct(std::vector<size_t> sz, Box<dim,T> dom)
	{
#ifdef DEBUG
		//! The size is wrong signal it

		if (sz.size() != dim)
		{std::cerr << "Error this factory has been specialized for catesian grid of dimension " << dim << "\n";}

#endif

		// Calculate the size of the hyper-cubes on each dimension

		T szd[dim];

		for (int i = 0 ; i < dim ; i++)
		{szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid<dim,void> g(sz);

		// Create a graph with the number of vertices equal to the number of
		// grid point

		//! Graph to construct

		Graph gp(g.size());

		/******************
		 *
		 * Create the edges and fill spatial
		 * information properties
		 *
		 ******************/

		//! Construct a key iterator

		grid_key_dx_iterator<dim> k_it(g);

		//! Iterate through all the elements

		while (k_it.isNext())
		{
			grid_key_dx<dim> key = k_it.get();

			// Vertex object
			auto obj = gp.vertex(g.LinId(key));

			// vertex spatial properties functor
			fill_prop<dim,T,decltype(gp.vertex(g.LinId(key))), typename to_boost_mpl<pos...>::type, sizeof...(pos) == 0 > flp(obj,szd,key);

			// fill properties

			boost::mpl::for_each< boost::mpl::range_c<int,0,sizeof...(pos)> >(flp);

			// Get the combinations of dimension d

			for (int d = dim-1 ; d >= dim_c ; d--)
			{
				// create the edges for that dimension

				std::vector<comb<dim>> c = hc.getCombinations_R(d);

				// for each combination calculate a safe linearization and create an edge

				for (int j = 0 ; j < c.size() ; j++)
				{
					// Calculate the element size

					T ele_sz = 0;

					// for each dimension multiply and reduce

					for (int s = 0 ; s < dim ; s++)
					{
						ele_sz += szd[s] * abs(c[j][s]);
					}

					// Calculate the end point vertex id
					// Calculate the start point id

					size_t start_v = g.LinId(key);
					size_t end_v = g.template LinId<CheckExistence>(key,c[j].getComb());

					// Add an edge and set the the edge property to the size of the face (communication weight)
					gp.template addEdge<CheckExistence>(start_v,end_v);
				}
			}

			// Fill vertex properties

			++k_it;
		}

		return gp;
	}
};


/*! \brief This class work as a functor
 *
 * For each number in the boost::mpl::vector (for example 3 6) set the properties of the vertex at the
 * specified id (3 6) with pos[d] * spacing[d] with d running from 0 to 1, pos[d] the position id of the vertex
 * spacing the grid spacing
 *
 * Example
 *
 * if we give a grid_key of dimension 2 4x4 the expression "pos[d] * spacing[d]"
 * will assume the value
 *
 * (0.0 0.0) (0.25 0.0) ...... (1.0 0.0)
 * (0.0 0.25)................. (1.0 0.25)
 * ....................................
 * (0.0 1.0).................. (1.0 1.0)
 *
 * and the properties 3 6 will be filled with the numbers 0.0 0.0    .......  1.0 1.0
 * progressively
 *
 * \tparam dim Dimensionality of the cartesian grid
 * \tparam dT type of the domain
 * \tparam G_v vertex type object
 * \tparam v boost::mpl::vector containing all the index to fill
 *
 */

template<unsigned int dim, typename dT, typename G_v, typename v>
class fill_prop<dim,dT,G_v,v,true>
{

public:

	//! Fill the object from where to take the properties
	fill_prop(G_v & g_v , const dT (& szd)[dim], grid_key_dx<dim> & gk)
	{}

	//! It call the function for each property we want to copy
    template<typename T>
    void operator()(T& t) const
    {}
};

/*! \brief This class construct a cartesian graph
 *
 * This class construct a cartesian graph
 *
 * \param dim dimensionality of the cartesian grid
 *
 */

template<unsigned int dim, typename Graph>
class CartesianGraphFactory
{

public:

	/*!
	 *
	 * \brief Construct a cartesian graph, with V and E edge properties
	 *
	 * Construct a cartesian graph, with V and E edge properties
	 *
	 * Each vertex is a subspace (Hyper-cube) of dimension dim, each vertex is
	 * connected with an edge if two vertex (Hyper-cube) share a element of dimension grater than
	 * dim_c
	 *
	 * \param sz Vector that store the size of the grid on each dimension
	 * \param dom Box enclosing the physical domain
	 *
	 * \tparam se Indicate which properties fill with the element weight. The
	 *           element weight is the point, line , surface, d dimensional object
	 *           in contact (in common between two hyper-cube). NO_EDGE indicate
	 *           no property will store this information
	 * \tparam T type of the domain like (int real complex ... )
	 * \tparam dim_c Connectivity dimension
	 * \tparam Memory class that create new memory
	 * \tparam pos... one or more integer indicating the spatial properties
	 *
	 */
	template <unsigned int se,typename T, unsigned int dim_c, int... pos>
	static Graph construct(std::vector<size_t> sz, Box<dim,T> dom)
	{
		return Graph_constructor_impl<dim,Graph,se,T,dim_c,pos...>::construct(sz,dom);
	}
};

#endif /* CARTESIANGRAPHFACTORY_HPP_ */
