/*
 * CartesianGraphFactory.hpp
 *
 *  Created on: Nov 28, 2014
 *      Author: i-bird
 */

#ifndef CARTESIANGRAPHFACTORY_HPP_
#define CARTESIANGRAPHFACTORY_HPP_

#include "Vector/map_vector.hpp"
#include "Graph/map_graph.hpp"
#include "Grid/grid_sm.hpp"
#include "Space/Shape/Box.hpp"
#include "Space/Shape/HyperCube.hpp"

#define NO_VERTEX_ID -1

/*! \brief Operator to fill the property 'prp' with the linearization of indexes
 *
 *  \tparam dim Dimension of the space
 *  \tparam G_v Graph
 *  \tparam prp Property to fill
 */
template<unsigned int dim, typename G_v, int prp>
struct fill_id
{
	//! function that fill with linearization indexes
	static inline void fill(G_v & g_v, const grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs)
	{
		g_v.template get<prp>() = gs.LinId(gk);
	}
};
/*! \brief Operator to fill the property in case there are no properties
 *
 *  \tparam dim Dimension of the space
 *  \tparam G_v Graph
 */
template<unsigned int dim, typename G_v>
struct fill_id<dim, G_v, NO_VERTEX_ID>
{
	//! function that fill with linearization indexes
	static inline void fill(G_v & g_v, const grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs)
	{
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
 * \tparam is_stub when is true, produce a trivial operator(),
 *         to use when v is an empty vector to avoid compilation error
 *
 */

template<unsigned int dim, int lin_id, typename dT, typename G_v, typename v, int impl>
class fill_prop
{
	//! Domain
	const Box<dim,dT> & domain;

	//! Reference to an array containing the spacing
	const dT (&szd)[dim];

	//! grid_key_dx Reference containing the actual position
	grid_key_dx<dim> & gk;

	//! Vertex object to fill
	G_v & g_v;

	//! grid info
	const grid_sm<dim, void> & gs;

public:

	//! Fill the object from where to take the properties
	fill_prop(G_v & g_v, const dT (&szd)[dim], grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs, const Box<dim,dT> & domain)
	:domain(domain), szd(szd), gk(gk), g_v(g_v), gs(gs)
	{
	}

	//! It call the function for each property we want to copy
	template<typename T>
	void operator()(T& t) const
	{
		typedef typename boost::fusion::result_of::at<v, boost::mpl::int_<T::value>>::type t_val;

		g_v.template get<t_val::value>() = gk.get(T::value) * szd[T::value] + domain.getLow(T::value);
		fill_id<dim, G_v, lin_id>::fill(g_v, gk, gs);
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

template<unsigned int dim, int lin_id, typename dT, typename G_v, typename v>
class fill_prop<dim, lin_id, dT, G_v, v, 0>
{

public:

	//! Fill the object from where to take the properties
	fill_prop(G_v & g_v, const dT (&szd)[dim], grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs, const Box<dim,dT> & domain)
	{
	}

	//! It call the function for each property we want to copy
	template<typename T>
	void operator()(T& t) const
	{
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

template<unsigned int dim, int lin_id, typename dT, typename G_v, typename v>
class fill_prop<dim, lin_id, dT, G_v, v, 2>
{
	//! Domain
	const Box<dim,dT> & domain;

	//! Reference to an array containing the spacing
	const dT (&szd)[dim];

	//! grid_key_dx Reference containing the actual position
	grid_key_dx<dim> & gk;

	//! Vertex object to fill
	G_v & g_v;

	//! grid info
	const grid_sm<dim, void> & gs;

public:

	//! Fill the object from where to take the properties
	fill_prop(G_v & g_v, const dT (&szd)[dim], grid_key_dx<dim> & gk, const grid_sm<dim, void> & gs, const Box<dim,dT> & domain)
	:domain(domain), szd(szd), gk(gk), g_v(g_v), gs(gs)
	{
	}

	//! It call the function for each property we want to copy
	template<typename T>
	void operator()(T& t) const
	{
		typedef typename boost::fusion::result_of::at<v, boost::mpl::int_<0>>::type t_val;
		typedef typename boost::mpl::at<typename G_v::T_type::type,t_val>::type s_type;

		for (size_t i = 0 ; i < std::extent<s_type>::value ; i++)
			g_v.template get<t_val::value>()[i] = 0.0;

		for (size_t i = 0 ; i < dim ; i++)
			g_v.template get<t_val::value>()[i] = gk.get(i) * static_cast<float>(szd[i]) + domain.getLow(i);

		fill_id<dim, G_v, lin_id>::fill(g_v, gk, gs);
	}
};

/*! \brief Operator for vector and scalar property
 *
 * \tparam i Size of the property
 * \tparam p Type of the property boost mpl
 * \tparam Graph Graph
 * \tparam pos Array of properties
 */
template<int i, typename p, typename Graph, int ... pos>
struct fill_prop_by_type
{

	//! Get the element 0
	typedef typename boost::mpl::at<p, boost::mpl::int_<0>>::type v_element;

	//! Get the property v_element (v_element is a number)
	typedef typename boost::mpl::at<typename Graph::V_type::type, v_element>::type pos_prop_type;

	enum
	{
		value = ((sizeof...(pos) != 0) * (std::is_array<pos_prop_type>::value + 1))
	};

};

/*! \brief Operator for vector and scalar property in the case there are no properties
 *
 * \tparam i Size of the property
 * \tparam p Type of the property
 * \tparam Graph Graph
 * \tparam pos Array of properties
 */
template<typename p, typename Graph, int ... pos>
struct fill_prop_by_type<0, p, Graph, pos...>
{
	enum
	{
		value = 0
	};

};

/*! \brief Graph constructor function specialization
 *
 * On C++ partial function specialization is not allowed, so we need a class to do it
 *
 * \see CartesianGraphFactory method construct
 *
 */

template<unsigned int dim, int lin_id, typename Graph, int se, typename T, unsigned int dim_c, int ... pos>
class Graph_constructor_impl
{
public:

	/*! \brief Construct a cartesian graph
	 *
	 * \param sz size of the partesian graph
	 * \param dom domain where this cartesian graph is defined (used to fill the coordinates)
	 * \param bc boundary conditions (torus or cube)
	 *
	 * \return the constructed graph
	 *
	 */
	static Graph construct(const size_t (& sz)[dim], Box<dim,T> dom, const size_t(& bc)[dim])
	{
		// Calculate the size of the hyper-cubes on each dimension
		T szd[dim];

		for (size_t i = 0; i < dim; i++)
		{
			szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];
		}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid_sm<dim, void> g(sz);

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

			typedef typename to_boost_vmpl<pos...>::type p;

			// vertex spatial properties functor

			fill_prop<dim, lin_id, T, decltype(gp.vertex(g.LinId(key))), typename to_boost_vmpl<pos...>::type, fill_prop_by_type<sizeof...(pos), p, Graph, pos...>::value> flp(obj, szd, key, g, dom);

			// fill properties

			boost::mpl::for_each<boost::mpl::range_c<int, 0, sizeof...(pos)> >(flp);

			// Get the combinations of dimension d

			for (long int d = dim-1 ; d >= dim_c ; d--)
			{
				// create the edges for that dimension

				std::vector<comb<dim>> c = hc.getCombinations_R(d);

				// for each combination calculate a safe linearization and create an edge

				for (size_t j = 0; j < c.size(); j++)
				{
					// Calculate the element size

					T ele_sz = 0;

					// for each dimension multiply and reduce


					for (size_t s = 0 ; s < dim ; s++)
						ele_sz += szd[s] * abs(c[j][s]);

					// Calculate the end point vertex id
					// Calculate the start point id

					size_t start_v = g.LinId(key);

					size_t end_v = g.template LinId<CheckExistence>(key,c[j].getComb(),bc);

					// Add an edge and set the the edge property to the size of the face (communication weight)
					gp.template addEdge<CheckExistence>(start_v, end_v).template get<se>() = ele_sz;
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

template<unsigned int dim, int lin_id, typename Graph, typename T, unsigned int dim_c, int ... pos>
class Graph_constructor_impl<dim, lin_id, Graph, NO_EDGE, T, dim_c, pos...>
{
public:

	/*! \brief Construct a cartesian graph
	 *
	 * \param sz size of the partesian graph
	 * \param dom domain where this cartesian graph is defined (used to fill the coordinates)
	 * \param bc boundary conditions (torus or cube)
	 *
	 * \return the constructed graph
	 *
	 */
	static Graph construct(const size_t ( & sz)[dim], Box<dim,T> dom, const size_t(& bc)[dim])
	{
		// Calculate the size of the hyper-cubes on each dimension

		T szd[dim];

		for (size_t i = 0; i < dim; i++)
		{
			szd[i] = (dom.getHigh(i) - dom.getLow(i)) / sz[i];
		}

		//! Construct an hyper-cube of dimension dim

		HyperCube<dim> hc;

		// Construct a grid info

		grid_sm<dim, void> g(sz);

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

			typedef typename to_boost_vmpl<pos...>::type p;

			// vertex spatial properties functor

			fill_prop<dim, lin_id, T, decltype(gp.vertex(g.LinId(key))), typename to_boost_vmpl<pos...>::type, fill_prop_by_type<sizeof...(pos), p, Graph, pos...>::value> flp(obj, szd, key, g, dom);

			// fill properties

			boost::mpl::for_each_ref<boost::mpl::range_c<int, 0, sizeof...(pos)> >(flp);

			// Get the combinations of dimension d

			for (long int d = dim-1 ; d >= dim_c ; d--)
			{
				// create the edges for that dimension

				std::vector<comb<dim>> c = hc.getCombinations_R(d);

				// for each combination calculate a safe linearization and create an edge

				for (size_t j = 0; j < c.size(); j++)
				{
					// Calculate the end point vertex id
					// Calculate the start point id

					size_t start_v = g.LinId(key);

					size_t end_v = g.template LinId<CheckExistence>(key,c[j].getComb(),bc);

					// Add an edge and set the the edge property to the size of the face (communication weight)
					gp.template addEdge<CheckExistence>(start_v, end_v);
				}
			}

			// Fill vertex properties

			++k_it;
		}

		return gp;
	}
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
	 * dim_c. One property can be used to store the contact size or the d-dimensional
	 * surface in common between two connected hyper-cube.
	 *
	 * \tparam se Indicate which properties fill with the contact size. The
	 *           contact size is the point, line , surface, d-dimensional object size
	 *           in contact (in common) between two hyper-cube. NO_EDGE indicate
	 *           no property will store this information
	 * \tparam id_prp property 'id' that stores the vertex id (with -1 it skip)
	 * \tparam T type of the domain like (int real complex ... )
	 * \tparam dim_c Connectivity dimension
	 * \tparam pos... (optional)one or more integer indicating the spatial properties
	 *
	 * \param sz store the size of the cartesian grid on each dimension
	 * \param dom Box enclosing the physical domain
	 * \param bc boundary conditions {PERIODIC = torus and NON_PERIODIC = cube}
	 *
	 * \return the constructed graph
	 *
	 */
	template<int se, int id_prp, typename T, unsigned int dim_c, int ... pos>
	static Graph construct(const size_t (&sz)[dim], Box<dim, T> dom, const size_t (& bc)[dim])
	{
		return Graph_constructor_impl<dim, id_prp, Graph, se, T, dim_c, pos...>::construct(sz, dom, bc);
	}
};

#endif /* CARTESIANGRAPHFACTORY_HPP_ */
